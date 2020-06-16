import numpy as np

N = 10  # Number of data points
D = 5  # Dimensionality, equivalent to the number of nodes
ETA = 0.2  # Sparsity: the probability of an edge being present
SEED = 0


def generate_graph(num_nodes, num_data, sparsity):
    """
    Create graph
    populate graph using sparsity
    ensure lower diagonal is empty
    """
    edges = D*(D - 1)/2
    np.random.seed(SEED)
    entries = np.random.choice([0, 1], size=int(edges), p=[ETA, 1-ETA])

    graph = np.zeros((num_nodes, num_nodes))

    edge_indexes = np.triu_indices(num_nodes, 1)
    graph[edge_indexes] = entries

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
    graph[np.where(graph)] = weight_values
    return graph

def simulate_data(weights, num_data, num_nodes):
    """
    for n in range(num_data):
        sample x_1 from a normal
    """
    data = np.zeros([num_data, num_nodes])
    
    for n in range(num_data):
        np.random.seed(SEED+n)
        noise = np.random.normal(0, 1, num_nodes)
        print(noise)
        data[n, 0] = noise[0]
        for i in range(1, num_nodes):
            data[n, i] = np.dot(data[n, :i], weights[:i, i]) + noise[i]
        print(data[n,:])
    print(data)
    return np.around(data, 3)


def main():
    graph = generate_graph(D, N, ETA)
    weights = create_weights(graph)
    data = simulate_data(weights, N, D)

    # Save file
    np.savetxt(
        "sim_data_seed_{}_N_{}_D_{}.csv".format(SEED, N, D),
        data,
        delimiter=','
    )

if __name__ == '__main__':
    main()
    



