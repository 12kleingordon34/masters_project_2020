import numpy as np

N = 1000  # Number of data points
D = 10  # Dimensionality, equivalent to the number of nodes
ETA = 0.3  # Sparsity: the probability of an edge being present
SEED = 0


def generate_graph(num_nodes, num_data, sparsity):
    """
    Create graph
    populate graph using sparsity
    ensure lower diagonal is empty
    """
    edges = D*(D - 1)/2
    np.random.seed(SEED)
    entries = np.random.choice([0, 1], size=int(edges), p=[1-ETA, ETA])
    graph = np.zeros((num_nodes, num_nodes))

    edge_indexes = np.triu_indices(num_nodes, 1)
    graph[edge_indexes] = entries
    print(graph)

    return graph

def create_weights(graph):
    """
    Number of all non-null nodes
    Sample a sequence of weights uniformally 
    Fill Graph with the weights
    """
    num_edges = np.sum(graph)
    np.random.seed(SEED)
    weight_values = np.random.uniform(-1, 1, int(num_edges))
    weights = np.zeros(graph.shape)
    weights[np.where(graph)] = weight_values
    return weights

def simulate_data(weights, num_data, num_nodes):
    """
    for n in range(num_data):
        sample x_1 from a normal
    """
    data = np.zeros([num_data, num_nodes])
    
    for n in range(num_data):
        np.random.seed(SEED+n)
        noise = np.random.normal(0, 1, num_nodes)
        data[n, 0] = noise[0]
        for i in range(1, num_nodes):
            data[n, i] = np.dot(data[n, :i], weights[:i, i]) + noise[i]
    return np.around(data, 3)


def main():
    graph = generate_graph(D, N, ETA)
    weights = create_weights(graph)
    data = simulate_data(weights, N, D)

    # Save file
    np.savetxt(
        "data/simulated_data/true_graph_seed_{}_N_{}_D_{}.csv".format(SEED, N, D),
        np.around(graph, 0),
        delimiter=','
    )
    np.savetxt(
        "data/simulated_data/true_weights_seed_{}_N_{}_D_{}.csv".format(SEED, N, D),
        np.around(weights, 3),
        delimiter=','
    )
    np.savetxt(
        "data/simulated_data/sim_data_seed_{}_N_{}_D_{}.csv".format(SEED, N, D),
        np.around(data, 3),
        delimiter=','
    )

if __name__ == '__main__':
    main()
