from __future__ import annotations

from backend.app.services.dataset_preprocessor import stratified_sample_cells


def test_stratified_sample_keeps_all_when_under_limit() -> None:
    group_values = ["A"] * 200 + ["B"] * 100
    cluster_values = ["c1"] * 150 + ["c2"] * 150

    selected, report = stratified_sample_cells(
        group_values=group_values,
        cluster_values=cluster_values,
        max_cells=500,
        seed=42,
    )

    assert len(selected) == 300
    assert report["mode"] == "full"
    assert report["input_cells"] == 300
    assert report["output_cells"] == 300


def test_stratified_sample_is_deterministic_and_bounded() -> None:
    group_values = ["A"] * 600 + ["B"] * 200
    cluster_values = (["a1"] * 500 + ["a2"] * 100) + (["b1"] * 150 + ["b2"] * 50)

    selected_1, report_1 = stratified_sample_cells(
        group_values=group_values,
        cluster_values=cluster_values,
        max_cells=500,
        seed=7,
    )
    selected_2, report_2 = stratified_sample_cells(
        group_values=group_values,
        cluster_values=cluster_values,
        max_cells=500,
        seed=7,
    )

    assert selected_1 == selected_2
    assert len(selected_1) == 500
    assert report_1["mode"] == "stratified_sampling"
    assert report_1["output_cells"] == 500
    assert report_1["group_counts_after"] == {"A": 375, "B": 125}
    assert report_2["group_counts_after"] == {"A": 375, "B": 125}
