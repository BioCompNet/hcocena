"""Community detection and scoring utilities."""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Tuple

from .data import ExpressionDataset
from .network import Graph


def detect_communities(graph: Graph, max_iterations: int = 25) -> List[List[str]]:
    """Detect communities using a deterministic label propagation algorithm."""

    labels = {node: idx for idx, node in enumerate(sorted(graph.nodes))}
    for _ in range(max_iterations):
        changes = 0
        for node in sorted(graph.nodes):
            neighbors = graph.neighbors(node)
            if not neighbors:
                continue
            neighbor_labels = [labels[nbr] for nbr in neighbors]
            counts = Counter(neighbor_labels)
            # choose the label with the highest count and break ties by label id
            best_label = min((-count, label) for label, count in counts.items())[1]
            if labels[node] != best_label:
                labels[node] = best_label
                changes += 1
        if changes == 0:
            break

    communities: Dict[int, List[str]] = defaultdict(list)
    for gene, label in labels.items():
        communities[label].append(gene)
    return [sorted(members) for members in communities.values()]


def communities_to_assignment(communities: Iterable[Iterable[str]]) -> Dict[str, int]:
    mapping: Dict[str, int] = {}
    for idx, community in enumerate(communities):
        for gene in community:
            mapping[gene] = idx
    return mapping


def compute_cluster_scores(dataset: ExpressionDataset, assignments: Dict[str, int]) -> List[Tuple[int, List[float]]]:
    """Average the expression values for each cluster."""

    grouped: Dict[int, List[List[float]]] = defaultdict(list)
    for gene, cluster in assignments.items():
        if gene in dataset.matrix:
            grouped[cluster].append(dataset.matrix[gene])

    scores: List[Tuple[int, List[float]]] = []
    for cluster, matrices in grouped.items():
        if not matrices:
            continue
        totals = [0.0] * len(dataset.samples)
        for vector in matrices:
            totals = [a + b for a, b in zip(totals, vector)]
        averages = [value / len(matrices) for value in totals]
        scores.append((cluster, averages))
    scores.sort(key=lambda item: item[0])
    return scores
