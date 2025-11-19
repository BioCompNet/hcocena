"""Graph helpers implemented without third-party dependencies."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, Iterable, List, Tuple

from .data import ExpressionDataset


@dataclass
class Graph:
    nodes: Dict[str, Dict[str, float]] = field(default_factory=dict)
    edges: Dict[Tuple[str, str], Dict[str, float]] = field(default_factory=dict)

    def add_node(self, name: str, **attributes) -> None:
        node = self.nodes.setdefault(name, {})
        node.update(attributes)

    def add_edge(self, first: str, second: str, **attributes) -> None:
        if first == second:
            return
        key = tuple(sorted((first, second)))
        edge = self.edges.setdefault(key, {})
        edge.update(attributes)
        self.add_node(first)
        self.add_node(second)

    def neighbors(self, node: str) -> List[str]:
        adjacent = []
        for (first, second) in self.edges:
            if first == node:
                adjacent.append(second)
            elif second == node:
                adjacent.append(first)
        return adjacent

    def degree(self, node: str) -> int:
        return len(self.neighbors(node))


def _pearson(x: List[float], y: List[float]) -> float:
    n = len(x)
    if n == 0:
        return 0.0
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    cov = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    var_x = sum((a - mean_x) ** 2 for a in x)
    var_y = sum((b - mean_y) ** 2 for b in y)
    if var_x == 0 or var_y == 0:
        return 0.0
    return cov / (var_x ** 0.5 * var_y ** 0.5)


def _spearman(x: List[float], y: List[float]) -> float:
    def rank(values: List[float]) -> List[float]:
        sorted_values = sorted(enumerate(values), key=lambda item: item[1])
        ranks = [0.0] * len(values)
        i = 0
        while i < len(sorted_values):
            j = i
            total_rank = 0.0
            while j < len(sorted_values) and sorted_values[j][1] == sorted_values[i][1]:
                total_rank += j + 1
                j += 1
            avg_rank = total_rank / (j - i)
            for k in range(i, j):
                ranks[sorted_values[k][0]] = avg_rank
            i = j
        return ranks

    rank_x = rank(x)
    rank_y = rank(y)
    return _pearson(rank_x, rank_y)


def _correlation(x: List[float], y: List[float], method: str) -> float:
    if method == "spearman":
        return _spearman(x, y)
    return _pearson(x, y)


def build_correlation_network(
    dataset: ExpressionDataset,
    method: str = "spearman",
    min_edge_weight: float = 0.5,
) -> Graph:
    graph = Graph()
    genes = list(dataset.matrix.keys())
    for gene in genes:
        values = dataset.matrix[gene]
        graph.add_node(
            gene,
            mean=sum(values) / len(values),
            variance=sum((value - sum(values) / len(values)) ** 2 for value in values) / len(values),
        )

    for gene_a, gene_b in combinations(genes, 2):
        weight = _correlation(dataset.matrix[gene_a], dataset.matrix[gene_b], method)
        if abs(weight) >= min_edge_weight:
            graph.add_edge(gene_a, gene_b, weight=weight)

    return graph


def integrate_networks(graphs: Iterable[Graph], strategy: str = "mean") -> Graph:
    graphs = list(graphs)
    if not graphs:
        raise ValueError("No graphs provided")

    integrated = Graph()
    accumulator: Dict[Tuple[str, str], List[float]] = {}

    for graph in graphs:
        for node, attrs in graph.nodes.items():
            integrated.add_node(node, **attrs)
        for edge, attrs in graph.edges.items():
            accumulator.setdefault(edge, []).append(attrs.get("weight", 0.0))

    for edge, weights in accumulator.items():
        if strategy == "sum":
            value = sum(weights)
        else:
            value = sum(weights) / len(weights)
        integrated.add_edge(edge[0], edge[1], weight=value, observations=len(weights))

    return integrated
