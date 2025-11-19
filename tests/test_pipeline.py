from hcocena_py import HCoCena
from hcocena_py.data import ExpressionDataset, zscore_normalize


def build_dataset() -> ExpressionDataset:
    samples = ["sample_a", "sample_b", "sample_c"]
    matrix = {
        "gene1": [5.0, 4.5, 5.5],
        "gene2": [2.0, 2.5, 1.5],
        "gene3": [3.0, 3.5, 3.5],
        "gene4": [4.0, 4.0, 3.5],
    }
    return ExpressionDataset(name="test", samples=samples, matrix=matrix)


def test_normalisation_zero_mean():
    dataset = build_dataset()
    normalised = zscore_normalize(dataset)
    for values in normalised.matrix.values():
        assert abs(sum(values)) < 1e-6


def test_workflow(tmp_path):
    workflow = HCoCena(min_edge_weight=0.2)
    workflow.add_dataset(build_dataset())
    workflow.build_all_networks()
    workflow.integrate_networks()
    workflow.detect_communities()
    scores = workflow.score_clusters()
    assert scores
    workflow.plot_cluster_heatmap(tmp_path / "heatmap.tsv")
    assert (tmp_path / "heatmap.tsv").exists()
    workflow.export_cluster_assignments(tmp_path / "clusters.csv")
    assert (tmp_path / "clusters.csv").exists()
