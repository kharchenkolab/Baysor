"""cellAdmix dataset plumbing for tabular (contract-format) molecule tables.

The stock Python bindings only build *Xenium* input stores
(``celladmix.dataset.CellAdmix.ensure_store`` raises for other formats), while
the C++ core has a full tabular store builder. The wheel we install (see
``install.sh``) exposes it as ``celladmix._core.build_tabular_store`` through a
small committed patch; :class:`TabularCellAdmix` wires that binding into the
normal ``CellAdmix`` lifecycle so ``fit()`` and ``fit.audit_admixture()`` work
unchanged on benchmark molecule tables.

``quick_cluster`` wraps the same C++ clustering the R API exposes as
``ds$cluster()`` — cellAdmix's quick cell typing used when no annotation is
supplied.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import celladmix as ca
from celladmix import _core


class TabularCellAdmix(ca.CellAdmix):
    """``CellAdmix`` subclass that builds tabular input stores.

    ``source`` is the molecules parquet/CSV path; ``schema`` maps the store
    columns (``x``, ``y``, ``gene``, ``cell``) onto the file's column names.
    """

    def __init__(self, molecules, *, output_dir, schema: dict | None = None, **kwargs):
        super().__init__(molecules, output_dir=output_dir, format="tabular", **kwargs)
        self.schema = {"x": "x", "y": "y", "gene": "gene", "cell": "cell"} | dict(schema or {})

    def ensure_store(self, *, force: bool = False,
                     materialize_molecules: bool = True, verbose: bool = True):
        if self.format != "tabular":
            return super().ensure_store(
                force=force, materialize_molecules=materialize_molecules, verbose=verbose)
        return _core.build_tabular_store(
            str(self.source),
            str(self.store_dir),
            x_col=self.schema["x"],
            y_col=self.schema["y"],
            gene_col=self.schema["gene"],
            cell_id_col=self.schema["cell"],
            keep_unassigned=self.keep_unassigned,
            materialize_molecules=materialize_molecules,
            force=force,
            num_threads=self.num_threads,
        )


def build_tabular_dataset(molecules_path: Path | str, output_dir: Path | str, *,
                          annotation=None, num_threads: int = 6,
                          schema: dict | None = None) -> TabularCellAdmix:
    """Create a tabular :class:`TabularCellAdmix` dataset object."""
    return TabularCellAdmix(
        molecules_path,
        output_dir=output_dir,
        schema=schema,
        annotation=annotation,
        num_threads=num_threads,
    )


def read_store_cells(store_dir: Path | str) -> list[str]:
    """Cell ids present in an input store."""
    cells = pd.read_parquet(Path(store_dir) / "cells.parquet", columns=["cell_id"])
    return cells["cell_id"].astype(str).tolist()


def quick_cluster(
    store_dir: Path | str,
    *,
    resolution: float = 1.0,
    seed: int = 1,
    min_molecules: int = 1,
    min_genes: int = 1,
    n_variable_genes: int = 1000,
    pca_dims: int = 30,
    graph_k: int = 15,
) -> pd.Series:
    """cellAdmix quick clustering as a cell-id-indexed type series.

    ``min_molecules``/``min_genes`` default to 1 and ``cells_max`` is
    unlimited so every store cell receives a label — the admixture audit
    treats cells without a type as a spurious pseudo-type.

    Clustering is forced to a single thread: the parallel graph step is
    racy (empirically 1 in 6 multi-threaded runs on the same store and seed
    produced a different partition, moving the total admixture rate by
    ~0.001), while single-threaded runs reproduce the partition exactly.
    The clustering itself takes well under a second on benchmark crops.
    """
    res = _core.cluster_store_counts(
        str(store_dir),
        min_molecules=int(min_molecules),
        min_genes=int(min_genes),
        cells_max=-1,
        n_variable_genes=int(n_variable_genes),
        pca_dims=int(pca_dims),
        graph_k=int(graph_k),
        cluster_resolution=float(resolution),
        compute_umap=False,
        umap_neighbors=15,
        umap_epochs=200,
        num_threads=1,
        umap_parallel_optimization=True,
        normalization_scale=5000.0,
        seed=int(seed),
    )
    labels = [f"cluster_{int(c)}" for c in res["cluster"]]
    return pd.Series(
        labels,
        index=np.asarray(res["cell_id"], dtype=str),
        name="celltype",
    )
