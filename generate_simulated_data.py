import numpy as np
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

N = 5000  # Number of data points
D = 4 # Dimensionality, equivalent to the number of nodes
ETA = 0.4  # Sparsity: the probability of an edge being present
SEED = 0
DISTRIBUTION = 'gaussian'


def generate_graph(num_nodes, num_data, sparsity, seed=SEED):
    """
    Create graph
    populate graph using sparsity
    ensure lower diagonal is empty
    """
    edges = D*(D - 1)/2
    np.random.seed(seed)
    entries = np.random.choice([0, 1], size=int(edges), p=[1-ETA, ETA])
    graph = np.zeros((num_nodes, num_nodes))

    edge_indexes = np.triu_indices(num_nodes, 1)
    graph[edge_indexes] = entries
    print(graph)

    return graph

def create_weights(graph, seed=SEED):
    """
    Number of all non-null nodes
    Sample a sequence of weights uniformally 
    Fill Graph with the weights
    """
    num_edges = np.sum(graph)
    np.random.seed(seed)
    weight_values = np.random.uniform(0.1, 2, int(num_edges))
    weights = np.zeros(graph.shape)
    weights[np.where(graph)] = weight_values
    return weights

def simulate_data(weights, num_data, num_nodes, distribution, seed):
    """
    for n in range(num_data):
        sample x_1 from a normal
    """
    data = np.zeros([num_data, num_nodes])
    
    for n in range(num_data):
        if distribution == 'gaussian':
            np.random.seed(seed+n)
            noise = np.random.normal(0, 1, num_nodes)
            data[n, 0] = noise[0] + 1
            for i in range(1, num_nodes):
                data[n, i] = np.dot(data[n, :], weights[:, i]) + noise[i]
        if distribution == 'gamma':
            data[n, 0] = np.random.gamma(shape=2, scale=2)
            for i in range(1, num_nodes):
                if np.dot(data[n, :i], weights[:i, i]) != 0:
                    data[n, i] = np.random.gamma(shape=2, scale=1/(np.dot(data[n, :i], weights[:i, i])))
                else:
                    data[n, i] = np.random.gamma(shape=2, scale=2)
    return np.around(data, 3)


def main():
    #graph = generate_graph(D, N, ETA, SEED)
    graph = np.array([
        [0,0,0,1,1,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,1,1,1,0]
    ])
    print(graph)
    weights = create_weights(graph, SEED)
    weights = graph * 2
    print(weights)
    print(graph.shape[0])
    data = simulate_data(weights, N, graph.shape[0], DISTRIBUTION, SEED)

    # Save file
    np.savetxt(
        "data/simulated_data/{}_true_graph_seed_{}_N_{}_D_{}.csv".format(DISTRIBUTION, SEED, N, D),
        graph.astype(int),
        delimiter=','
    )
    np.savetxt(
        "data/simulated_data/{}_true_weights_seed_{}_N_{}_D_{}.csv".format(DISTRIBUTION, SEED, N, D),
        np.around(weights, 3),
        delimiter=','
    )
    np.savetxt(
        "data/simulated_data/{}_sim_data_seed_{}_N_{}_D_{}.csv".format(DISTRIBUTION, SEED, N, D),
        np.around(data, 3),
        delimiter=','
    )
    np.savetxt(
        "data/simulated_data/bootstrap_sample_N_{}.csv".format(N),
        np.around(data, 3),
        delimiter=','
    )

#if __name__ == '__main__':
#    main()
