"""
Code Sample for MBTA Data Analyst Co-Op role
Jamie Koo

"""

import networkx as nx 
import pandas as pd
import matplotlib.pyplot as plt

class Graph:
    def __init__(self, vertices=None, edges=None):
        ''' Initialize an empty graph with vertices and edges '''
        self.graph = {v: set() for v in (vertices or [])}
        for u, v in (edges or []):
            self.add_edge(u, v)

    def add_vertex(self, vertex):
        ''' Add a vertex to the graph '''
        self.graph.setdefault(vertex, set())

    def add_edge(self, u, v):
        ''' Add an edge to the graph and ensure both vertices exist '''
        self.add_vertex(u)
        self.add_vertex(v)
        self.graph[u].add(v)
        self.graph[v].add(u)

    def remove_edge(self, u, v):
        ''' Remove an edge from the graph '''
        if v in self.graph.get(u, set()):
            self.graph[u].remove(v)
            self.graph[v].remove(u)

    def __getitem__(self, vertex):
        ''' Get neighbors of a vertex '''
        return self.graph.get(vertex, set())

    def __repr__(self):
        ''' Get a string representation of the graph '''
        return '\n'.join([f'[{v}] => {neighbors}' for v, neighbors in self.graph.items()])

    def is_adjacent(self, u, v):
        ''' Check if two vertices are adjacent '''
        return v in self.graph.get(u, set())

    def get_vertices(self):
        ''' Get a list of all vertices in the graph '''
        return list(self.graph.keys())

    def get_edges(self):
        ''' Get a list of all edges in the graph '''
        return [(u, v) for u in self.graph for v in self.graph[u] if u < v]

    def num_vertices(self):
        ''' Get the number of vertices in the graph '''
        return len(self.graph)

    def num_edges(self):
        ''' Get the number of edges in the graph '''
        return sum(len(neighbors) for neighbors in self.graph.values()) // 2 

    def deg(self, vertex):
        ''' Get the degree of a vertex '''
        return len(self.graph.get(vertex, set()))

    def degree_distribution(self):
        ''' Calculate the degree distribution of the graph '''
        distribution = {}
        for v in self.graph:
            degree = self.deg(v)
            distribution[degree] = distribution.get(degree, 0) + 1
        return distribution

    def to_df(self):
        ''' Convert graph to a Pandas DataFrame containing the edges '''
        edges = self.get_edges()
        return pd.DataFrame(edges, columns=['u', 'v'])

    def from_df(self, df):
        ''' Load a graph from a Pandas DataFrame containing unique vertices
        and edges '''
        vertices = pd.unique(df['u'].append(df['v']))
        edges = list(df.itertuples(index=False, name=None))
        self.__init__(vertices, edges)

    def visualize(self, fig=1, directed=False, bipartite_nodes = []):
        ''' Plot the graph using networkx and matplotlib '''
        df = self.to_df()
        options = {
            "font_size": 10,
            "node_size": 500,
            "node_color": "white",
            "edgecolors": "black",
            "linewidths": 1,
            "width": 1,
        }
        graph_type = nx.Graph()
        G = nx.from_pandas_edgelist(df, df.columns[0], df.columns[1], create_using=graph_type)
        if len(bipartite_nodes) > 0:
            pos = nx.bipartite_layout(G, bipartite_nodes)
        else:
            pos = nx.circular_layout(G)
        plt.figure(figsize = (10, 10))
        nx.draw_networkx(G, pos, with_labels=True, **options)
        plt.show()
        

class WeightedGraph(Graph):
    def __init__(self, vertices=None, edges=None):
        ''' Initialize a weighted graph with vertices and weighted edges '''
        super().__init__(vertices)
        self.weights = {}
        for u, v, w in (edges or []):
            self.add_edge(u, v, w)

    def add_edge(self, u, v, weight=None):
        ''' Add a weighted edge to the graph '''
        super().add_edge(u, v)
        self.weights.setdefault(u, {})[v] = weight
        self.weights.setdefault(v, {})[u] = weight

    def __getitem__(self, vertex):
        ''' Get all vertices adjacent to a vertex with weights '''
        return self.weights.get(vertex, {})

    def to_df(self):
        ''' Convert the graph to a Pandas Dataframe with edges and weights '''
        edges = [(u, v, w) for u in self.weights for v, w in self.weights[u].items()]
        return pd.DataFrame(edges, columns=['u', 'v', 'w'])

    def from_df(self, df):
        ''' Load graph from Pandas Dataframe with vertices and weight '''
        super().from_df(df[['u', 'v']])
        for _, row in df.iterrows():
            self.add_edge(row['u'], row['v'], row['w'])

    def visualize(self, fig = 1, directed=False):
        '''' Render the graph using networkx and matplotlib libraries '''
        df = self.to_df()
        options = {
            "font_size": 24,
            "node_size": 1000,
            "node_color": "white",
            "edgecolors": "black",
            "linewidths": 3,
            "width": 3,
        }
        graph_type = nx.Graph()
        G = nx.from_pandas_edgelist(df, df.columns[0], df.columns[1],
                                    edge_attr=True, create_using=graph_type)
        pos = nx.spring_layout(G)
        plt.figure(fig)
        nx.draw_networkx(G, pos, with_labels=True, **options)
        labels = nx.get_edge_attributes(G,'w')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
        plt.show()
        
    def subgraph(self, vertices):
        ''' Create a subgraph containing only specified vertices '''
        edges = [(u, v, w) for u in vertices for v, w in self[u].items()]
        return WeightedGraph(vertices, edges)

    def from_csv(self, filename):
        ''' Load a graph from a CSV file '''
        df = pd.read_csv(filename, usecols=["gene", "disease", "num_positive"])
        for _, row in df[df['num_positive'] > 0].iterrows():
            self.add_edge(row['gene'], row['disease'], int(row['num_positive']))
    
    def asthma_analysis(self):
        ''' Perform analysis of genes linked to asthma '''
        asthma_genes = [v for v, neighbors in self.graph.items() if 'asthma' in neighbors]
        asthma_linked_diseases = {gene: list(self.graph[gene]) for gene in asthma_genes}
        subg = Graph()
        for gene, diseases in asthma_linked_diseases.items():
            subg.add_vertex(gene)
            for disease in diseases:
                subg.add_edge(gene, disease)
        return subg, asthma_genes
    
def plot_distribution(degree_distribution):
    ''' Plot the degree distribution of a graph on a log scale '''
    degrees, counts = zip(*degree_distribution.items())
    plt.scatter(degrees, counts)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('GAD degree distribution')
    plt.xlabel('Degree')
    plt.ylabel('Number of Nodes')
    plt.grid(True)
    plt.show()
        
def main():
    # Part 1: Testing Graph class
    V = list("ABCDEFGH")
    E = [('A', 'B'), ('A', 'C'), ('A', 'G'), ('A', 'H'),
         ('B', 'C'), ('B', 'F'), ('C', 'D'), ('D', 'E'),
         ('E', 'F'), ('H', 'F')]
    g = Graph(V, E)
    print('\nUndirected Graph')
    print('|V|:', g.num_vertices())
    print('|E|:', g.num_edges())
    print('Adjacent to A:', g['A'])
    print(g)

    # Part 2: Testing WeightedGraph class
    V = list("ABCDEFGH")
    E = [('A', 'B', 5), ('A', 'C', 3), ('A', 'G', 2), ('A', 'H', 9),
         ('B', 'C', 0), ('B', 'F', 3), ('C', 'D', 1),
         ('D', 'E', 12), ('E', 'F', 16), ('H', 'F', 8)]
    g = WeightedGraph(V, E)
    print('\nUndirected Weighted Graph')
    print('|V|:', g.num_vertices())
    print('|E|:', g.num_edges())
    print('Adjacent to A:', g['A'])
    print(g)
    subg = g.subgraph(['A', 'B', 'C'])
    print(subg)
    
    # Part 3: Actual Code - GAD Data Analysis
    wg = WeightedGraph()
    wg.from_csv('gad_data.csv')
    degree_distribution = wg.degree_distribution()
    plot_distribution(degree_distribution)
    print(degree_distribution)

    # Asthma Analysis
    print("\nPerforming Asthma Analysis...")
    subg, asthma_gene = wg.asthma_analysis()
    print('Asthma-linked subgraph:')
    print(subg)
    subg.visualize(bipartite_nodes = asthma_gene)

if __name__ == '__main__':
    main()
