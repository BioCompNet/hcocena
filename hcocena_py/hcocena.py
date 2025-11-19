"""Stateful Python workflow that mirrors the R package."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from .clustering import compute_cluster_scores, communities_to_assignment, detect_communities
from .data import ExpressionDataset, load_expression_table, zscore_normalize
from .network import Graph, build_correlation_network, integrate_networks
from .visualization import save_heatmap_table


@dataclass
class HCoCena:
    min_edge_weight: float = 0.5
    correlation_method: str = "spearman"
    datasets: Dict[str, ExpressionDataset] = field(default_factory=dict)
    networks: Dict[str, Graph] = field(default_factory=dict)
    integrated_network: Optional[Graph] = None
    _assignments: Optional[Dict[str, int]] = None
    _cluster_scores: Optional[List[Dict[str, float]]] = None

    def add_dataset(self, dataset: ExpressionDataset, normalise: bool = True) -> None:
        if dataset.name in self.datasets:
            raise ValueError(f"Dataset '{dataset.name}' already exists")
        if normalise:
            dataset = zscore_normalize(dataset)
        self.datasets[dataset.name] = dataset

    def add_dataset_from_file(self, path: str | Path, name: Optional[str] = None) -> None:
        dataset = load_expression_table(path, name)
        self.add_dataset(dataset)

    def build_network(self, dataset_name: str) -> Graph:
        dataset = self.datasets[dataset_name]
        graph = build_correlation_network(
            dataset,
            method=self.correlation_method,
            min_edge_weight=self.min_edge_weight,
        )
        self.networks[dataset_name] = graph
        return graph

    def build_all_networks(self) -> Dict[str, Graph]:
        return {name: self.build_network(name) for name in self.datasets}

    def integrate_networks(self) -> Graph:
        if not self.networks:
            self.build_all_networks()
        self.integrated_network = integrate_networks(self.networks.values())
        return self.integrated_network

    def detect_communities(self) -> Dict[str, int]:
        if self.integrated_network is None:
            self.integrate_networks()
        communities = detect_communities(self.integrated_network)
        self._assignments = communities_to_assignment(communities)
        return self._assignments

    def score_clusters(self) -> List[Dict[str, float]]:
        if self._assignments is None:
            self.detect_communities()

        results: List[Dict[str, float]] = []
        for dataset in self.datasets.values():
            scores = compute_cluster_scores(dataset, self._assignments)
            for cluster_id, values in scores:
                row: Dict[str, float] = {"dataset": dataset.name, "cluster": cluster_id}
                for sample, value in zip(dataset.samples, values):
                    row[sample] = value
                results.append(row)
        self._cluster_scores = results
        return results

    def get_cluster_assignments(self) -> Dict[str, int]:
        if self._assignments is None:
            self.detect_communities()
        return dict(self._assignments)

    def export_cluster_assignments(self, path: str | Path) -> None:
        assignments = self.get_cluster_assignments()
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            handle.write("gene,cluster\n")
            for gene, cluster in sorted(assignments.items()):
                handle.write(f"{gene},{cluster}\n")

    def plot_cluster_heatmap(self, output_path: str | Path) -> None:
        scores = self._cluster_scores if self._cluster_scores is not None else self.score_clusters()
        summary: Dict[int, Dict[str, float]] = {}
        for row in scores:
            cluster = int(row["cluster"])
            dataset = str(row["dataset"])
            values = [value for key, value in row.items() if key not in {"cluster", "dataset"}]
            if not values:
                continue
            summary.setdefault(cluster, {})[dataset] = sum(values) / len(values)
        save_heatmap_table(summary, output_path)

    def summary(self) -> List[Dict[str, int]]:
        assignments = self.get_cluster_assignments()
        rows: List[Dict[str, int]] = []
        for dataset in self.datasets.values():
            covered = len(set(dataset.matrix) & set(assignments))
            rows.append(
                {
                    "dataset": dataset.name,
                    "genes": len(dataset.matrix),
                    "samples": len(dataset.samples),
                    "covered_genes": covered,
                }
            )
        return rows
