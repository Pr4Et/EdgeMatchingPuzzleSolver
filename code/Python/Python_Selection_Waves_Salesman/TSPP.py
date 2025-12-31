import networkx as nx
import random
import matplotlib.pyplot as plt


# Parameters
N = 4
M = 6
DELAY_RANGE = (1, 10)

# Initialize directed graph
G = nx.DiGraph()

# Generate M random edges with random delays
edges = set()
while len(edges) < M:
    u = random.randint(0, N - 1)
    v = random.randint(0, N - 1)
    if u != v:
        edge_pair = tuple(sorted((u, v)))  # ensures (u,v) and (v,u) are treated the same
        edges.add(edge_pair)

# Add edges in both directions with the same delay
for u, v in edges:
    delay = random.randint(*DELAY_RANGE)
    G.add_edge(u, v, delay=delay)
    G.add_edge(v, u, delay=delay)

# Messages: { (origin, destination): [origin_node] }
messages = {}
counters = {}
for u, v in G.edges():
    messages[(u, v)] = [u]
    counters[(u, v)] = G[u][v]['delay']

arrival_lists = {node: [] for node in G.nodes()}

# Simulation loop
while True:
    # Step 1: Reset arrival lists and decrement counters
    for node in arrival_lists:
        arrival_lists[node] = []

    for msg_id in list(counters.keys()):
        counters[msg_id] -= 1
        if counters[msg_id] == 0:
            arrival_lists[msg_id[1]].append(msg_id)

    # Step 2: Process arrivals
    for pnode in G.nodes():
        arrivals = arrival_lists[pnode]

        # Handle arrival lists with more than one message
        if len(arrivals) > 1:
            chosen = random.choice(arrivals)
            arrivals.remove(chosen)
            for msg_id in arrivals:
                counters[msg_id] += 1

        # After checking for >1, now check for exactly 1 message
        if len(arrivals) == 1:
            chosen = arrivals[0]
            messages[chosen].append(pnode)
            temp=messages[chosen]
            # Loop detection and termination condition
            if messages[chosen].count(chosen[1]) > 1:
                if set(messages[chosen]).issuperset(set(G.nodes())):  #to compare with the Ising model expression, any closure with the path should do
                    print("Final message path:", messages[chosen])

                    path_nodes = messages[chosen]  # e.g., [1, 3, 5]
                    path_edges = list(zip(path_nodes, path_nodes[1:]))  # Result: [(1, 3), (3, 5)]
                    # Compute layout
                    pos = nx.spring_layout(G)  # Or use your preferred layout

                    # Draw base graph
                    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray')

                    # Highlight path edges
                    nx.draw_networkx_edges(G, pos, edgelist=path_edges, edge_color='red', width=2)

                    # Annotate each edge in path with order number
                    for idx, edge in enumerate(path_edges, start=1):
                        # Midpoint between nodes
                        x = (pos[edge[0]][0] + pos[edge[1]][0]) / 2
                        y = (pos[edge[0]][1] + pos[edge[1]][1]) / 2
                        plt.text(x, y, str(idx), fontsize=12, color='black', ha='center', va='center')

                    plt.title('Graph Visualization with Traversal Edge Order')
                    plt.show()

                    exit()
                else:
                    del messages[chosen]
                    del counters[chosen]
                    continue

            # Propagate to connected nodes
            origin = pnode
            prev_origin = chosen[0]
            for neighbor in G.successors(origin):
                if neighbor != prev_origin:
                    new_key = (origin, neighbor)
                    messages[new_key] = messages[chosen].copy()
                    counters[new_key] = G[origin][neighbor]['delay']

            # Delete the chosen message after propagation
            del messages[chosen]
            del counters[chosen]