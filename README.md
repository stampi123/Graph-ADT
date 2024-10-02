# Graph-ADT

## Overview
This repository contains an implementation of a weighted, undirected graph in C++. The graph supports various operations, such as finding connected components, calculating shortest paths (both weighted and unweighted), and computing the smallest connecting threshold between two nodes. It also includes a Disjoint Set data structure for union-find operations.

## Features
- **Disjoint Set (Union-Find)**: Efficiently manages node sets and supports union and find operations.
- **Graph Construction**: Loads graph data from a CSV file containing edge lists (node pairs with weights).
- **Basic Operations**:
  - Count the number of nodes and edges.
  - Get the neighbors of a node.
  - Calculate the weight of an edge between two nodes.
- **Graph Algorithms**:
  - Unweighted shortest path using Breadth-First Search (BFS).
  - Weighted shortest path using Dijkstra's algorithm.
  - Identify connected components based on a given edge weight threshold.
  - Compute the smallest connecting threshold between two nodes.
