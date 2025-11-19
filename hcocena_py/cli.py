"""Command line interface for running the Python workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from .hcocena import HCoCena


def parse_dataset_argument(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Dataset specification must use the NAME=PATH format")
    name, path = value.split("=", 1)
    return name, Path(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the Python implementation of hCoCena on expression matrices")
    parser.add_argument(
        "--dataset",
        action="append",
        type=parse_dataset_argument,
        required=True,
        help="Dataset specification in the form NAME=PATH",
    )
    parser.add_argument("--output", type=Path, default=Path("hcocena_output"), help="Output directory")
    parser.add_argument("--min-edge", type=float, default=0.5, help="Minimum absolute correlation for edges")
    parser.add_argument("--method", type=str, default="spearman", help="Correlation method")
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    workflow = HCoCena(min_edge_weight=args.min_edge, correlation_method=args.method)

    for name, path in args.dataset:
        workflow.add_dataset_from_file(path, name)

    workflow.build_all_networks()
    workflow.integrate_networks()
    workflow.detect_communities()
    scores = workflow.score_clusters()

    output_dir = args.output
    output_dir.mkdir(parents=True, exist_ok=True)

    with (output_dir / "cluster_scores.tsv").open("w", encoding="utf-8") as handle:
        if scores:
            header = [key for key in scores[0].keys()]
            handle.write("\t".join(map(str, header)) + "\n")
            for row in scores:
                handle.write("\t".join(str(row.get(column, "")) for column in header) + "\n")

    assignments = workflow.get_cluster_assignments()
    workflow.export_cluster_assignments(output_dir / "gene_clusters.csv")
    workflow.plot_cluster_heatmap(output_dir / "cluster_heatmap.tsv")


if __name__ == "__main__":  # pragma: no cover
    main()
