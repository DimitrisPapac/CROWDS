#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from simulate import *


def main():
    phi = 0.75                                            # Probability to forward.
    graph_file = "sample_graph.txt"                       # The path to the file containing the graph.
    corrupted_file = "./sample_corrupted_users_file.txt"  # The path to the file containing the IDs of the corrupt users.
    users_file = "./sample_users_file.txt"                # The path to the file containing the initiators for each simulation.
    broken_paths = 5                                      # Maximum number of paths that the adversary is allowed to destroy.
    fix_strategy = "initiator"                            # The fix-strategy employed when a path is destroyed (either 'initiator' or 'last-honest').

    corrupted_users = parse_corrupt_users_file(corrupted_file)
    initiators = parse_initiators_file(users_file)

    for initiator in initiators:
        if initiator in corrupted_users:  # skip corrupt users
            continue

        graph = parse_graph_file(graph_file)
        for i in range(len(graph)):
            graph[i][i] = 1
        result = run(phi, initiator, graph, corrupted_users, broken_paths, fix_strategy)
        print(f"User {result} delivers message to web server.")


if __name__ == '__main__':
    main()
