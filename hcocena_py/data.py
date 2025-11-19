"""Data handling utilities for the pure Python hCoCena workflow."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional


def _infer_separator(path: Path) -> str:
    if path.suffix in {".tsv", ".txt"}:
        return "\t"
    return ","


@dataclass
class ExpressionDataset:
    """Stores a gene expression matrix."""

    name: str
    samples: List[str]
    matrix: Dict[str, List[float]]
    metadata: Optional[Dict[str, str]] = None

    def copy(self) -> "ExpressionDataset":
        return ExpressionDataset(
            name=self.name,
            samples=list(self.samples),
            matrix={gene: list(values) for gene, values in self.matrix.items()},
            metadata=None if self.metadata is None else dict(self.metadata),
        )


def load_expression_table(path: str | Path, name: Optional[str] = None) -> ExpressionDataset:
    """Load an expression table from ``path``.

    The first column must contain gene identifiers and subsequent columns the
    sample measurements.
    """

    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(file_path)

    separator = _infer_separator(file_path)
    with file_path.open("r", newline="") as handle:
        reader = csv.reader(handle, delimiter=separator)
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise ValueError("Expression file must contain at least one sample column")
        samples = header[1:]
        matrix: Dict[str, List[float]] = {}
        for row in reader:
            if len(row) != len(header):
                raise ValueError("Row length does not match header length")
            gene = row[0]
            if gene in matrix:
                raise ValueError(f"Duplicate gene identifier: {gene}")
            values = [float(value) for value in row[1:]]
            matrix[gene] = values

    dataset_name = name or file_path.stem
    return ExpressionDataset(name=dataset_name, samples=samples, matrix=matrix)


def zscore_normalize(dataset: ExpressionDataset) -> ExpressionDataset:
    """Return a z-score normalised copy of the dataset."""

    from math import sqrt

    normalised = {}
    for gene, values in dataset.matrix.items():
        mean = sum(values) / len(values)
        variance = sum((value - mean) ** 2 for value in values)
        std = sqrt(variance / len(values)) if variance else 0.0
        if std == 0:
            normalised[gene] = [0.0 for _ in values]
        else:
            normalised[gene] = [(value - mean) / std for value in values]

    return ExpressionDataset(name=dataset.name, samples=list(dataset.samples), matrix=normalised, metadata=dataset.metadata)


def iterate_matrix(matrix: Dict[str, List[float]]) -> Iterable[tuple[str, List[float]]]:
    for gene, values in matrix.items():
        yield gene, list(values)
