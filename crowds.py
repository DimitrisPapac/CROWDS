#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import Union, List, Set

ALLOWABLE_STRATEGIES = {"last-honest", "initiator"}


# A function for parsing a file containing the graph's adjacency matrix.
def parse_graph_file(filepath: Union[str, bytes]) -> List[List[int]]:
    """
    A function for parsing a file containing the network's graph.

    Parameters:
        filepath (Union[str, bytes]) -- the path to the file containing the network's graph.

    Returns:
        List[List[int]] -- the adjacency matrix representing the network of jondos.
    """
    with open(filepath, mode="r", encoding="utf-8") as f:
        rows = f.readlines()

        graph = []
        num_of_columns = []
        max_num_of_columns = 0
        for row in rows:
            tokens = [int(x) for x in row.split()]
            num_of_columns.append(len(tokens))
            if num_of_columns[-1] > max_num_of_columns:
                max_num_of_columns = num_of_columns[-1]

            graph.append([int(x) for x in tokens])
            # Note: Generator object allows for short-circuit evaluation of any.
            if any((x != 0 and x != 1) for x in graph[-1]):
                raise Exception("Graph file entries must be binary (i.e., 0 or 1).")

        for i in range(len(graph)):
            graph[i] += (max_num_of_columns - num_of_columns[i]) * [0]   # pad zeros if necessary

        if len(graph) != len(graph[0]):
            raise Exception("Adjacency matrix must be square.")

    return graph


def parse_corrupt_users_file(filepath: Union[str, bytes]) -> Set[int]:
    """
    A function for parsing a file containing the corrupt users' identities.

    Parameters:
        filepath (Union[str, bytes]) -- the path to the file containing the corrupt users' identities.

    Returns:
        Set[int] -- The set of the ids of all corrupt jondos in the network.
    """
    with open(filepath, mode="r", encoding="utf-8") as f:
        rows = f.readlines()

        jondo_ids = set()   # ignores duplicates
        for row in rows:
            tokens = [x for x in row.split()]
            if len(tokens) != 1:
                raise Exception("Corrupt users file must contain one entry per line")
            jondo_id = int(tokens[0])
            jondo_ids.add(jondo_id)

    return jondo_ids


def parse_initiators_file(filepath: Union[str, bytes]) -> List[int]:
    """
    A function for parsing a file containing the initiators' identities for each simulation.

    Parameters:
        filepath (Union[str, bytes]) -- the path to the file containing the initiator's identities identities.

    Returns:
        List[int] -- A list coprising of the jondo ids that will serve as initiators in each protocol run.
    """
    with open(filepath, mode="r", encoding="utf-8") as f:
        rows = f.readlines()

        initiators = []
        for row in rows:
            tokens = [x for x in row.split()]
            if len(tokens) != 1:
                raise Exception("Only one initiator per simulation. Check your <users-file>.")
            initiator = int(tokens[0])
            initiators.append(initiator)

    return initiators


def is_symmetric(a: List[List[int]]) -> bool:
    """
    A function for checking if a matrix is symmetric.

    Parameters:
        a (List[List[int]]) -- the matrix in question, represented as a list of lists.

    Returns:
        bool -- True iff the input matrix is symmetric.
    """
    m, n = len(a), len(a[0])

    # Note: List comprehensions would prevent short-circuiting from taking place.
    return all(a[i][j] == a[j][i] for i in range(m) for j in range(i, n))


# A function for performing a single simulation run of CROWDS for the given parameters.
def run(phi: float,
        initiator: int,
        graph: List[List[int]],
        corrupted_users: Set[int],
        broken_paths: int = 0,
        fix_strategy:str = "last-honest") -> int:
    """
    A function for performing a single simulation run of CROWDS for the given parameters.

    Parameters:
        phi (float) -- The (constant) probability of forwarding a request to a neighboring jondo.
        initiator (int) -- The jondo from which the protocol starts.
        graph (List[List[int]]) -- The adjacency matrix representing the network.
        corrupted_users (Set[int]) -- The set of corrupted users in the network.

    Keyword arguments:
        broken_paths (int) -- The maximum allowable number of destroyed paths (default: 0).
        fix_strategy (str) -- The protocol's fixing strategy. (default: "last-honest").

    Returns:
        int -- The id of the jondo that delivers the message to the web server.
    """

    assert (0 < phi < 1)  # phi is a probability. Moreover, for phi in {0, 1} the simulation is meaningless.
    if phi <= 0.5:
        print(f"Warning: Forwarding parameter {phi = }. Probable innocence does not hold.")

    m = len(graph)            # number of vertices in the graph.
    c = len(corrupted_users)  # number of corrupt users (collaborators).
    assert (type(broken_paths) == int and 0 <= broken_paths <= m * (m - 1) / 2)

    if m < (c + 1) * (1 + 0.5 / (phi - 0.5)):
        print("Warning: Inequality m >= (c + 1) * (1 + 0.5/(phi - 0.5)) does not hold. \
        Probable innocence does not hold.")

    assert fix_strategy in ALLOWABLE_STRATEGIES

    from random import uniform, randint

    current_node = initiator  # initiator must not be corrupt
    while True:
        r = uniform(0, 1)
        if r < phi:  # forward message to a random neighboring jondo (possibly itself)
            row = graph[current_node]
            neighbors = [j for j in range(len(row)) if row[j] == 1]   # indices of connected nodes
            next_node = neighbors[randint(0, len(neighbors) - 1)]     # or nextNode = random.choice(neighbors)

            # Note: Adversary cannot prevent a jondo from forwarding a message to itself.
            if next_node != current_node and next_node in corrupted_users and broken_paths > 0:
                broken_paths -= 1  # decrease counter
                graph[current_node][next_node] = graph[next_node][current_node] = 0  # sever/destroy path
                if fix_strategy == "initiator":
                    current_node = initiator
                # else: # fix_strategy == "last-honest":
                #    currentNode = currentNode
            else:
                current_node = next_node
        else:   # deliver message to server
            return current_node


# TEST:
if __name__ == "__main__":  # If the module is run directly as a standalone program (i.e., not imported)
    import sys

    phi = float(sys.argv[1])         # Probability to forward.
    graph_file = sys.argv[2]         # The path to the file containing the graph.
    corrupted_file = sys.argv[3]     # The path to the file containing the IDs of the corrupt users.
    users_file = sys.argv[4]         # The path to the file containing the initiators for each simulation.
    broken_paths = int(sys.argv[5])  # Maximum number of paths that the adversary is allowed to destroy.
    fix_strategy = sys.argv[6]       # The fix-strategy employed when a path is destroyed.

    corrupted_users = parse_corrupt_users_file(corrupted_file)
    initiators = parse_initiators_file(users_file)

    for initiator in initiators:
        if initiator in corrupted_users:
            continue
        graph = parse_graph_file(graph_file)
        for i in range(len(graph)):
            graph[i][i] = 1
        result = run(phi, initiator, graph, corrupted_users, broken_paths, fix_strategy)
        print(f"User {result} delivers message to web server.")
