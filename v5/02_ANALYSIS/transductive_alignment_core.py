from __future__ import annotations

from typing import Any

import numpy as np


def _symmetric_matrix_power(matrix: np.ndarray, exponent: float) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    safe_eigenvalues = np.clip(eigenvalues, np.finfo(float).eps, None)
    return (eigenvectors * np.power(safe_eigenvalues, exponent)) @ eigenvectors.T


def coral_align_source_to_target(
    source: np.ndarray, target: np.ndarray, regularization: float = 1e-6
) -> np.ndarray:
    """Align source mean and covariance to an unlabeled target batch using CORAL."""
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    if source.ndim != 2 or target.ndim != 2 or source.shape[1] != target.shape[1]:
        raise ValueError("source and target must be two-dimensional with the same feature count")
    if source.shape[0] < 2 or target.shape[0] < 2:
        raise ValueError("source and target require at least two rows")
    if not np.isfinite(source).all() or not np.isfinite(target).all():
        raise ValueError("CORAL requires finite source and target values")
    if regularization <= 0:
        raise ValueError("regularization must be positive")

    source_mean = source.mean(axis=0, keepdims=True)
    target_mean = target.mean(axis=0, keepdims=True)
    centered_source = source - source_mean
    centered_target = target - target_mean
    feature_count = source.shape[1]
    source_covariance = centered_source.T @ centered_source / (source.shape[0] - 1)
    target_covariance = centered_target.T @ centered_target / (target.shape[0] - 1)
    ridge = regularization * np.eye(feature_count)
    transform = _symmetric_matrix_power(source_covariance + ridge, -0.5) @ _symmetric_matrix_power(
        target_covariance + ridge, 0.5
    )
    return centered_source @ transform + target_mean


def assert_label_free_alignment_contract(
    *,
    target_labels_used_for_alignment: bool,
    target_labels_used_for_model_selection: bool,
    target_labels_used_for_assessment: bool,
) -> dict[str, Any]:
    """Reject analyses that access target diagnoses before predictions are locked."""
    if target_labels_used_for_alignment or target_labels_used_for_model_selection:
        raise ValueError("target labels must not be used for alignment or model selection")
    if not target_labels_used_for_assessment:
        raise ValueError("target labels are required only for final locked assessment")
    return {
        "target_labels_used_for_alignment": False,
        "target_labels_used_for_model_selection": False,
        "target_labels_used_only_for_assessment": True,
        "deployment_scope": "batch_level_transductive_stress_test",
    }
