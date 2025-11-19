"""Text based summaries that mimic a heatmap-like overview."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List


def save_heatmap_table(matrix: Dict[int, Dict[str, float]], output_path: str | Path) -> None:
    """Persist the numeric summary to ``output_path`` in TSV format."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    datasets: List[str] = sorted({dataset for row in matrix.values() for dataset in row})
    with output.open("w", encoding="utf-8") as handle:
        handle.write("cluster\t" + "\t".join(datasets) + "\n")
        for cluster in sorted(matrix):
            values = [f"{matrix[cluster].get(dataset, 0.0):.4f}" for dataset in datasets]
            handle.write(f"{cluster}\t" + "\t".join(values) + "\n")
