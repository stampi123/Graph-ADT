#include "Graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <utility>
#include <queue>
#include <algorithm>
#include <set>
#include <limits>
#include <unordered_set>

using namespace std;

class DisjointSet {
private:
    unordered_map<string, string> parent;
    unordered_map<string, int> rank;

public:
    DisjointSet() {}

    void makeSet(const string& u) {
        parent[u] = u;
        rank[u] = 0;
    }

    string find(const string& u) {
        if (parent[u] != u)
            parent[u] = find(parent[u]);
        return parent[u];
    }

    void unionSets(const string& u, const string& v) {
        string rootU = find(u);
        string rootV = find(v);
        if (rootU == rootV) return;

        if (rank[rootU] < rank[rootV])
            parent[rootU] = rootV;
        else if (rank[rootU] > rank[rootV])
            parent[rootV] = rootU;
        else {
            parent[rootV] = rootU;
            rank[rootU]++;
        }
    }
};


struct Node {
    string label;
    vector<pair<string, double>> neighbors;  
};

unordered_map<string, Node> nodes_;

struct CompareDist {
    bool operator()(const pair<double, string>& a, const pair<double, string>& b) {
        return a.first > b.first;  
    }
};

Graph::Graph(const char* const & edgelist_csv_fn) {
 ifstream edge_file(edgelist_csv_fn);
    if (!edge_file.is_open()) {
        cerr << "Error: Unable to open edge list file '" << edgelist_csv_fn << "'." << endl;
        return;
    }

    string line;
    while (getline(edge_file, line)) {
        istringstream ss(line);
        string u_label, v_label, weight_str;
        if (getline(ss, u_label, ',') && getline(ss, v_label, ',') && getline(ss, weight_str, ',')) {
            double weight = stod(weight_str); 

           
            if (nodes_.count(u_label) == 0) {
                nodes_[u_label] = Node{u_label, {}};
            }
            nodes_[u_label].neighbors.push_back(make_pair(v_label, weight));

           
            if (nodes_.count(v_label) == 0) {
                nodes_[v_label] = Node{v_label, {}};
            }
            nodes_[v_label].neighbors.push_back(make_pair(u_label, weight));
        } else {
            cerr << "Error: Invalid edge format: " << line << endl;
        }
    }
}

unsigned int Graph::num_nodes() {
    // TODO
      std::set<std::string> unique_nodes;

  
  for (const auto& node_pair : nodes_) {
    unique_nodes.insert(node_pair.first);
  }

  
  return unique_nodes.size();
}

vector<string> Graph::nodes() {
    // TODO
  vector<string> all_nodes;
  for (auto const& node_pair : nodes_) {
    all_nodes.push_back(node_pair.first);
  }
  return all_nodes;
}

unsigned int Graph::num_edges() {
    // TODO
    unsigned int edge_count = 0;
    for (std::pair<const std::string, Node> const& node_pair : nodes_) {
        edge_count += node_pair.second.neighbors.size();
    }
    return edge_count/2;
}

unsigned int Graph::num_neighbors(string const & node_label) {
    // TODO
     if (nodes_.count(node_label) == 0) {
        return 0; 
    }
    return nodes_[node_label].neighbors.size();
}

double Graph::edge_weight(string const & u_label, string const & v_label) {
    // TODO
      if (nodes_.count(u_label) == 0 || nodes_.count(v_label) == 0) {
        return -1; 
    }

    for (const auto& neighbor : nodes_[u_label].neighbors) {
        if (neighbor.first == v_label) {
            return neighbor.second;
        }
    }

    
    return -1;
}

vector<string> Graph::neighbors(string const & node_label) {
    // TODO
  if (nodes_.count(node_label) == 0) {
    return vector<string>(); 
  }
  vector<string> neighbor_labels;
  for (auto const& neighbor : nodes_[node_label].neighbors) {
    neighbor_labels.push_back(neighbor.first);
  }
  return neighbor_labels;
}

vector<string> Graph::shortest_path_unweighted(string const & start_label, string const & end_label) {
    // TODO
      
  unordered_map<string, bool> visited;
  for (auto const& node_pair : nodes_) {
    visited[node_pair.first] = false;
  }

  
  queue<string> bfs_queue;

 
  unordered_map<string, string> parent;


  bfs_queue.push(start_label);
  visited[start_label] = true;

  while (!bfs_queue.empty()) {
    string current_node = bfs_queue.front();
    bfs_queue.pop();

    
    if (current_node == end_label) {
      break;
    }

    
    for (auto const& neighbor : neighbors(current_node)) {
      if (!visited[neighbor]) {

        bfs_queue.push(neighbor);

        visited[neighbor] = true;
        parent[neighbor] = current_node;
      }
    }
  }

  
  if (visited[end_label]) {
    vector<string> path;
    string current = end_label;


    while (current != start_label) {
      path.push_back(current);
      current = parent[current];
    }
    path.push_back(start_label);
    reverse(path.begin(), path.end()); 
    return path;
  }

  
  return vector<string>();

}

vector<tuple<string,string,double>> Graph::shortest_path_weighted(string const & start_label, string const & end_label) {
    // TODO
 priority_queue<pair<double, string>, vector<pair<double, string>>, CompareDist> pq;

    
    unordered_map<string, double> distance;
    unordered_map<string, string> parent;

    
    for (const auto& node_pair : nodes_) {
        distance[node_pair.first] = numeric_limits<double>::infinity();
    }

    
    distance[start_label] = 0;
    pq.push({0, start_label});

    // Dijkstra's algorithm
    while (!pq.empty()) {
        auto current = pq.top();
        pq.pop();
        double current_distance = current.first;
        string u = current.second;

        
        if (u == end_label) {
            break;
        }

        
        for (const auto& neighbor : neighbors(u)) {
            
            string v = neighbor;
            double weight = edge_weight(u, v);

            
            if (weight >= 0 && current_distance + weight < distance[v]) {
                
                distance[v] = current_distance + weight;
                parent[v] = u;
                pq.push({distance[v], v});
            }
        }
    }

    
    vector<tuple<string, string, double>> path;
    if (distance[end_label] < numeric_limits<double>::infinity()) {
        string current = end_label;
        while (current != start_label) {
            string prev = parent[current];
            double edge_weight = this->edge_weight(prev, current);
            if (edge_weight >= 0) {  
                path.push_back(make_tuple(prev, current, edge_weight));
            }
            current = prev;
        }
        reverse(path.begin(), path.end()); 
    }

    
    if (start_label == end_label) {
        path.push_back(make_tuple(start_label, end_label, -1));
    }

    return path;
}

vector<vector<string>> Graph::connected_components(double const & threshold) {
    // TODO
    unordered_map<string, bool> visited;
  vector<vector<string>> components;

  
  for (const auto& node_pair : nodes_) {
    string current_node = node_pair.first;

    
    if (visited[current_node]) {
      continue;
    }

  
    vector<string> component;

    
    queue<string> bfs_queue;
    bfs_queue.push(current_node);
    visited[current_node] = true;

    
    while (!bfs_queue.empty()) {
      current_node = bfs_queue.front();
      bfs_queue.pop();

      
      component.push_back(current_node);

     
      for (const auto& neighbor : neighbors(current_node)) {
        if (!visited[neighbor]) {
          
          double edge_weight = this->edge_weight(current_node, neighbor);
          if (edge_weight <= threshold) {
            bfs_queue.push(neighbor);
            visited[neighbor] = true;
          }
        }
      }
    }

 
    components.push_back(component);
  }

  return components;
}

double Graph::smallest_connecting_threshold(string const & start_label, string const & end_label) {
    // TODO
     if (start_label == end_label) {
        return 0; 
    }

    DisjointSet ds;

   
    for (const auto& node_pair : nodes_) {
        ds.makeSet(node_pair.first);
    }

    vector<tuple<string, string, double>> edges_sorted;
    for (const auto& node_pair : nodes_) {
        const string& u = node_pair.first;
        for (const auto& neighbor : node_pair.second.neighbors) {
            const string& v = neighbor.first;
            double weight = neighbor.second;
            edges_sorted.push_back(make_tuple(u, v, weight));
        }
    }
    sort(edges_sorted.begin(), edges_sorted.end(), [](const tuple<string, string, double>& a, const tuple<string, string, double>& b) {
        return get<2>(a) < get<2>(b);
    });

    
    double smallest_threshold = -1;

   
    for (const auto& edge : edges_sorted) {
        string u = get<0>(edge);

        string v = get<1>(edge);
        double weight = get<2>(edge);

        
        ds.unionSets(u, v);

       
        if (ds.find(start_label) == ds.find(end_label)) {

            smallest_threshold = weight;
            break; 
        }
    }

    return smallest_threshold;
}