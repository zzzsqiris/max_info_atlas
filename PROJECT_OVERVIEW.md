# max_info_atlas project overview and script audit

Audit date: 2026-08-08

## Reorganization update

The whole-embryo experiment folders were subsequently consolidated under
`MouseWholeEmbryo_MapInfo/`. The active shared reference is now in
`MouseWholeEmbryo_MapInfo/reference/whole_embryo/`; current spatial datasets
are under `MouseWholeEmbryo_MapInfo/datasets/`; and the rejected STARmap,
MERFISH EC-only, and MERFISH whole-embryo experiments are preserved under
`MouseWholeEmbryo_MapInfo/archive/`. The detailed script audit below describes
the pre-reorganization source folders and remains useful as provenance, but
those old top-level paths are no longer current entry points.

This document is the persistent orientation note for future Codex sessions. It
describes the repository as inspected in the working tree, including the
dataset folders that are currently untracked. It does not claim that those
untracked files are committed or part of the current branch history.

## Repository status and scope

The repository is `/Users/shiqizhong/Code/max_info_atlas`. The current branch is
`ovarian-xenium-reproduction`. The last committed change is the configurable,
deterministic Leiden seed. At audit time, the working tree already contained
user modifications to the ovarian plotting scripts and `src/`, plus untracked
dataset folders, configs, and `DOWNLOADS.md`. Those changes predate this note
and must be preserved.

The formal package is under `src/max_info_atlases/`. The dataset folders are
experiment and provenance layers around that package; they are not alternate
implementations of the core Map Information algorithm.

## Structure

```text
max_info_atlas/
├── src/max_info_atlases/          Core reusable package and max-info CLI
├── config/                        Core and ovarian run YAML files
├── XeniumOvarianReproduction/    Current active ovarian Xenium workflow
├── LabelTransfer/                 Previous STARmap + scRNA label-transfer pilot
├── LabelTransfer_E15.5/           E15.5 MERFISH endothelial-cell workflow
├── LabelTransfer_E15.5_WholeEmbryo/
│                                  Full E15.5 MERFISH / Qiu reference workflow
├── MouseWholeEmbryo_VisiumHD/    Visium HD conversion and retained data
├── MouseWholeEmbryo_StereoSeq/   Placeholder; no scripts or data at audit time
├── tests/                         Unit and small end-to-end/UGE test harnesses
├── README.md                      Package installation and CLI overview
├── TESTING.md                     Test and pipeline examples
├── DOWNLOADS.md                   Removed-data provenance inventory
├── setup.py, pyproject.toml       Packaging and Cython build configuration
└── max_info_env.yml               Conda environment specification
```

`data/` and `results/` under the dataset folders are generally ignored or
large working artifacts. Raw and processed data should not be committed. The
`.DS_Store` and `__pycache__` files seen during the audit are generated files,
not project scripts.

## Where the core pipeline lives

The reusable pipeline is in `src/max_info_atlases/` and is exposed through the
`max-info` console command configured in `setup.py`.

- `features/`: expression, PCA, local-frequency, derived-feature extraction,
  and k-nearest-neighbor graph construction.
- `clustering/`: Leiden, PhenoGraph, K-means, and LDA implementations.
- `percolation/`: edge-list management, graph percolation, scoring, and score
  aggregation.
- `run/`: YAML parsing, manifests, expected-output tracking, job generation,
  dependency-aware submission, status, and error analysis.
- `uge/`: Hoffman2/UGE job-list, chunk, submission, monitoring, and template
  helpers.
- `comparison/`: ARI comparison workflow.
- `cli/main.py`: interactive commands, workers, run orchestration, and result
  aggregation.

The core configured workflow is:

```text
AnnData + Sections.npy + section_XY.npy
    -> max-info run prepare
    -> features
    -> graphs
    -> clustering
    -> percolation
    -> aggregation / score reduction
    -> dataset-specific plots and interpretation
```

`RunOrchestrator.submit_all()` uses the ordered steps `features`, `graphs`,
`clustering`, `percolation`, and `aggregation`, with each later UGE step held
on the previous step. The normal entry points for a configured run are:

```bash
max-info run summary --config <config.yaml>
max-info run prepare --config <config.yaml>
max-info run submit --config <config.yaml> --step <step>
max-info run status --config <config.yaml>
```

For an ordinary expression-based run, the scientific defaults represented by
the current configs are PCA50, a correlation kNN graph with `k=15`, a 50-point
log-spaced Leiden resolution sweep, `max_k=500`, and `pbond_steps=101`.
Newer configs explicitly set Leiden `random_seed: 42`; older results without
that setting must not automatically be treated as identical deterministic
runs.

## Dataset-specific folders and current datasets

### `XeniumOvarianReproduction` — current active workflow

This is the active branch’s dataset workflow. It contains a 407,124-cell
ovarian Xenium gene-only AnnData input, exported `Sections.npy`, and full and
refined Map Information results. `config/ovarian_full_map_info.yaml` is the
full run; `config/ovarian_local_smoke.yaml` is a 10,000-cell local smoke test;
`config/ovarian_refined_map_info.yaml` and
`config/ovarian_refined_map_info_seed42.yaml` are refined sweeps that reuse
full-run features, graphs, and edge lists.

Typical dataset-specific commands are:

```bash
python XeniumOvarianReproduction/scripts/prepare_map_info_inputs.py
max-info run prepare --config config/ovarian_full_map_info.yaml
max-info run submit --config config/ovarian_full_map_info.yaml --step <step>
python XeniumOvarianReproduction/scripts/plot_resolution_sweep.py
python XeniumOvarianReproduction/scripts/plot_full_map_info_results.py
python XeniumOvarianReproduction/scripts/plot_voronoi_clusters.py
```

The preparation script only exports spatial inputs and an optional contiguous
smoke-test subset; it does not normalize expression or modify the source
AnnData. The ovarian plotting scripts are reporting tools, not core pipeline
workers.

### `MouseWholeEmbryo_VisiumHD`

This folder contains 10x Visium HD segmented mouse-embryo raw archives,
barcode-coordinate metadata, and a converted
`mouse_embryo_visium_hd_raw.h5ad`. `scripts/convert_visium_hd.py` is the
dataset-specific converter. There is currently no root run config or complete
Map Information workflow for this dataset, so the converter is provenance and
rebuild code rather than a current core-pipeline entry point.

### `MouseWholeEmbryo_StereoSeq`

The folder is currently empty. It is a dataset placeholder, not an implemented
pipeline.

### `LabelTransfer` — previous STARmap pilot

This folder contains conversion and inspection scripts for a STARmap BZ9
spatial sample and GSE116470 frontal-cortex scRNA reference, plus an earlier
one-shot correlation-kNN label-transfer script and pilot plots. The source
datasets and generated outputs were removed; see `DOWNLOADS.md` for accessions
and provenance. These scripts are reproducibility code for the previous pilot,
not current `src/` core code.

### `LabelTransfer_E15.5` — E15.5 EC branch

This is the smaller E15.5 MERFISH endothelial-cell workflow. It converts the
two scRNA replicates and MERFISH region, optionally selects EC candidates,
prepares shared-gene matrices and UGE chunks, computes correlation-kNN
neighbors, combines 50 reference Leiden resolutions by majority vote, exports
Map Information labels/spatial inputs, and optionally builds a scRNA-imputed
expression branch. `WORKFLOW.md` is the authoritative task-specific run order.

The task configs cover scRNA reference clustering, MERFISH ST-only, transferred
labels, and imputed-expression Map Information.

### `LabelTransfer_E15.5_WholeEmbryo` — full E15.5 branch

This is a later full-MERFISH/Qiu-reference branch. It converts the Qiu run-20
reference, clusters the 215,284-cell E15.5 reference, transfers labels to the
219,525-cell MERFISH target while preserving source-row mappings, and supports
full transfer, ST-only, and chunked scRNA-imputation Map Information runs.

Its scripts are largely a dataset-specific superset of the smaller E15.5
branch: it handles zero-variance target cells, full-target unassigned labels,
larger resources, and disk-safe imputation chunks. It has no committed
workflow document at audit time, so its shell wrappers and configs are the
available operational record.

## Standard workflows by stage

### 1. Convert or prepare the AnnData input

Use the converter matching the source format. Do not substitute converters
across datasets merely because they have similar names.

| Dataset/input | Typical script | Output role |
|---|---|---|
| STARmap + GSE116470 | `LabelTransfer/scripts/convert_to_anndata.py` | Previous pilot AnnData files |
| E15.5 scRNA + MERFISH region | `LabelTransfer_E15.5/scripts/convert_to_anndata.py` | E15.5 raw AnnData files |
| Qiu whole-embryo Matrix Market | `LabelTransfer_E15.5_WholeEmbryo/scripts/convert_qiu_e15_5_reference.py` | Filtered sparse reference AnnData |
| Visium HD segmented output | `MouseWholeEmbryo_VisiumHD/scripts/convert_visium_hd.py` | Visium HD raw AnnData |
| Ovarian Xenium | `XeniumOvarianReproduction/scripts/prepare_map_info_inputs.py` | Sections and XY arrays for existing AnnData |

Converters validate shape, names, sparsity, counts, coordinates, or metadata
alignment and often write through a temporary output. They are intentionally
dataset-specific because the input formats and biological filtering rules
differ.

### 2. Optional label transfer

For the E15.5 branches, the normal order is:

```text
prepare_label_transfer.py
    -> submit_label_transfer.py [preview first, then --submit]
    -> check_label_transfer.py
    -> combine_label_transfer.py
    -> prepare_mapinfo_inputs.py
```

The mathematics is shared-gene library-size normalization to 10,000, `log1p`,
exact correlation kNN with `k=15`, and majority voting for primary labels.
Weighted labels and vote/distance diagnostics are retained for comparison, not
used as the primary Map Information labels in the documented E15.5 workflow.

The full-embryo branch additionally uses
`target_source_indices.npy` and writes `unassigned_zero_shared_expression` for
target cells excluded from correlation kNN because their retained shared-gene
variance is zero.

### 3. Optional expression imputation

The E15.5 branch uses `impute_expression_from_neighbors.py`; the full branch
uses `impute_full_expression_from_neighbors.py`,
`check_imputation_chunks.py`, and
`merge_imputation_chunks_on_disk.py`. The full branch reuses the verified
label-transfer neighbor chunks, averages normalized full reference expression
with inverse-distance weights, validates every chunk, and then sends the
merged/imputed AnnData through the normal core feature/graph/clustering/
percolation/aggregation pipeline.

### 4. Run Map Information

Prepare the task-specific `Sections.npy`, XY arrays, preexisting label arrays,
and, where appropriate, the shared `maxK=500` spatial edge list. Then use the
core `max-info run` commands with the appropriate YAML config. Dataset shell
wrappers must not be mistaken for replacements for the core worker code.

### 5. Plot and interpret

Plotting scripts consume reduced score CSV files, cluster arrays, and spatial
coordinates after a run. They do not calculate the Map Information score.
For large Voronoi plots, preserve explicit coordinate/label length checks and
do not silently downsample a final result.

## Project conventions and assumptions

- AnnData row order is the primary alignment key. Spatial coordinates are
  expected in `adata.obsm["spatial"]` with shape `(n_cells, 2)`.
- Core per-section inputs use `Sections.npy` and
  `xy_coordinates/<section>_XY.npy`; the section filename stem must match the
  section label used by the pipeline.
- Raw count matrices remain sparse where possible. Avoid densifying large
  expression matrices except in the deliberately bounded label-transfer
  shared-gene step.
- Label-transfer features and labels are not interchangeable with core
  clustering features. Preserve the source-row mapping whenever a subset is
  created.
- Primary transferred labels are majority-vote labels. Weighted labels are
  diagnostics unless a task config explicitly says otherwise.
- Config paths are not uniformly portable. In particular, refined ovarian
  configs contain absolute Hoffman2 paths for reused feature/graph/edge-list
  directories and need review before use on another machine.
- `max-info run prepare` creates job lists, chunks, and a manifest; workers are
  `max-info run-features`, `run-graph`, `run-clustering`, `run-percolation`,
  and `run-aggregation`.
- Core tests should be run with project dependencies installed. If an unrelated
  pytest plugin prevents collection, use the established isolated invocation
  `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 pytest ...`.
- Do not infer a successful HPC run from a submitted job ID. Check output
  files, manifests, logs, score counts, and row/coordinate alignment.
- Raw/processed datasets and generated results are intentionally outside the
  normal Git-tracked source surface. Check Git status before making provenance
  claims.

## Audit of scripts outside `src/`

The categories below are primary recommendations, not proof that a script has
never been run. `Core` means needed by the reusable package/test or the current
ovarian workflow; `Dataset-specific` means needed only to reproduce that task;
`Duplicate` means substantially overlapping implementation; `Experimental /
obsolete` means an archive candidate after its results and provenance are
confirmed. A script can have more than one label.

### Root and tests

| Script | Classification and recommendation |
|---|---|
| `setup.py` | **Core** packaging/Cython build entry point; keep. |
| `tests/create_test_dataset.py` | **Core test support**; keep for reproducible small inputs. |
| `tests/run_end_to_end_test.sh` | **Core test harness** for local execution; keep. |
| `tests/run_uge_pipeline.sh` | **Core test harness** for UGE execution; keep if Hoffman2 regression testing remains supported. |
| `tests/test_clustering.py`, `tests/test_path_utils.py`, `tests/test_run_error_analysis.py` | **Core tests**; keep. |

### `LabelTransfer/scripts`

| Script | Classification and recommendation |
|---|---|
| `convert_to_anndata.py` | STARmap/GSE116470 **dataset-specific** converter; keep only for previous-pilot reproducibility. |
| `inspect_gse116470.py`, `inspect_starmap.py` | Read-only **experimental/provenance** checks; archive after recording the checks if the pilot is closed. |
| `transfer_labels_knn.py` | Earlier one-shot STARmap **dataset-specific** transfer; it is not the chunked 50-resolution E15.5 workflow. Archive candidate if STARmap is no longer an active comparison. |
| `plot_gse115746_umap.py`, `plot_starmap_mapinfo.py` | Previous-pilot **dataset-specific plotting**; archive candidates. `plot_starmap_mapinfo.py` reuses plotting helpers from the ovarian folder. |

### `LabelTransfer_E15.5/scripts`

| Script | Classification and recommendation |
|---|---|
| `convert_to_anndata.py` | E15.5 raw 10x/MERFISH **dataset-specific** converter; keep if the E15.5 branch is reproducible. |
| `prepare_merfish_ec_subset.py` | E15.5 biological filtering **dataset-specific** step; keep with the branch. |
| `prepare_label_transfer.py` | E15.5 chunked transfer preparation; **dataset-specific** and substantially duplicated by the WholeEmbryo superset. Keep until commonization. |
| `run_label_transfer_chunk.py` | Substantially identical to the WholeEmbryo copy; **Duplicate**. Use one shared copy, preferably based on the newer WholeEmbryo version. |
| `submit_label_transfer.py` | Dataset-specific UGE wrapper; **Duplicate** structure with WholeEmbryo, but resource/job-name defaults differ. Keep as a thin parameterized wrapper or commonize. |
| `check_label_transfer.py` | E15.5 completeness check; **Dataset-specific**, with a smaller feature set than the WholeEmbryo superset. Keep if this branch remains active. |
| `combine_label_transfer.py` | E15.5 all-resolution vote/summary writer; **Dataset-specific**, overlapping the WholeEmbryo copy. Keep for the branch until shared logic is extracted. |
| `prepare_mapinfo_inputs.py` | E15.5 transferred-label and spatial/edge-list preparation; **Dataset-specific**. Keep with the branch. |
| `prepare_scrna_clustering_sections.py` | One-off reference section-array preparation; **Dataset-specific** and simple. Archive after the reference artifact is permanently retained. |
| `impute_expression_from_neighbors.py` | E15.5 EC imputation branch; **Dataset-specific**. Keep only if the imputed comparison is still required. |
| `build_imputed_expression_job.sh` | One E15.5 UGE imputation wrapper; **Dataset-specific**. Keep with the imputation branch or archive after output/provenance capture. |
| `plot_label_transfer_confidence.py` | Transfer diagnostic plotting; **Dataset-specific**. Keep for confidence/distance review; otherwise archive. |
| `plot_mapinfo_sweep.py` | Older sweep plotting copy; **Duplicate** of the WholeEmbryo sweep implementation. Archive after adopting the newer parameterized copy. |
| `plot_mapinfo_voronoi.py` | Older Voronoi plotting copy; **Duplicate** of the WholeEmbryo superset. Archive after adopting the newer copy. |
| `plot_st_only_comparison.py` | Three-way E15.5 comparison plotting; **Dataset-specific/experimental**. Keep only while ST-only versus transfer versus imputation is an active analysis. |

### `LabelTransfer_E15.5_WholeEmbryo/scripts`

| Script | Classification and recommendation |
|---|---|
| `convert_qiu_e15_5_reference.py` | Qiu Matrix Market converter with strict filtering/validation; **Dataset-specific**, keep for full-embryo reproducibility. |
| `convert_qiu_e15_5_reference_job.sh` | One-off conversion UGE wrapper; **Experimental/provenance**. Archive after conversion is verified and source provenance is recorded. |
| `prepare_reference_clustering_sections.py` | Full-reference section-array preparation; **Dataset-specific**, keep if reference clustering must be rebuilt. |
| `build_reference_features_job.sh` | Manual reference feature job with shape checks; **Dataset-specific**. Keep while rebuilding the reference; archive as a one-off after the artifact is stable. |
| `prepare_label_transfer.py` | Full-embryo transfer preparation; **Dataset-specific and Duplicate**, a newer superset of the smaller E15.5 version. Use its zero-variance/source-index logic as the basis for a shared implementation. |
| `run_label_transfer_chunk.py` | Functionally identical to the smaller E15.5 copy; **Duplicate**. Preferred existing canonical copy for commonization. |
| `submit_label_transfer.py` | Full-embryo UGE wrapper; **Duplicate** structure with resource differences. Parameterize rather than maintain two copies. |
| `check_label_transfer.py` | Full-embryo completeness check with source-index validation; **Duplicate/superset**. Preferred basis for a shared checker. |
| `combine_label_transfer.py` | Full-target vote/summary writer with unassigned cells and transactional output; **Duplicate/superset**. Preferred basis for shared transfer-combine logic. |
| `prepare_mapinfo_inputs.py` | Full transferred-label Map Information input/edge-list preparation; **Dataset-specific**. Keep for the full-transfer branch. |
| `prepare_st_only_full_mapinfo.py` | Full-target ST-only spatial/edge-list preparation; **Dataset-specific** comparison branch. Keep only while ST-only is active. |
| `prepare_imputation_mapinfo_inputs.py` | Imputed-data spatial/section export; **Dataset-specific** and separate from transferred-label preparation. Keep with imputation. |
| `impute_full_expression_from_neighbors.py` | Full-embryo chunked imputation worker; **Dataset-specific**. Keep if the imputation branch is retained. |
| `check_imputation_chunks.py` | Required read-only chunk coverage/value validator; **Dataset-specific**, keep with imputation. |
| `merge_imputation_chunks_on_disk.py` | Required disk-safe imputation merge; **Dataset-specific**, keep with imputation. |
| `build_full_imputed_expression_array_job.sh` | Full imputation array wrapper; **Dataset-specific**. Keep if reruns on UGE are expected. |
| `build_full_imputed_expression_pilot_job.sh` | Leading-cell pilot wrapper; **Experimental**. Archive after pilot behavior/resources are documented. |
| `combine_label_transfer_job.sh` | Two-command UGE convenience wrapper; **Dataset-specific/duplicate orchestration**. Prefer explicit checker then combiner commands or a parameterized wrapper. |
| `prepare_label_transfer_job.sh` | Preparation plus inline invariant checks; **Dataset-specific**. Keep if it remains the documented reproducibility command. |
| `plot_mapinfo_sweep.py` | Newer, argument-driven sweep plotter; **Duplicate** of the smaller E15.5 copy and preferred existing canonical plotting copy. |
| `plot_mapinfo_voronoi.py` | Newer Voronoi plotter with dataset titles and source-index support; **Duplicate/superset** and preferred existing copy. |
| `run_full_mapinfo_pipeline.sh` | Convenience polling wrapper around the core CLI; **Experimental/obsolete candidate** because it relies on textual status matching and hardcoded paths. Verify outputs before archiving. |

### Other dataset folders

| Script | Classification and recommendation |
|---|---|
| `MouseWholeEmbryo_VisiumHD/scripts/convert_visium_hd.py` | Visium HD **dataset-specific** converter with sparse gene collapsing, coordinate validation, and transactional output; keep as rebuild/provenance code if Visium HD remains in scope. |
| `XeniumOvarianReproduction/scripts/prepare_map_info_inputs.py` | **Core for the current ovarian dataset**; keep as the current spatial-input preparation entry point. |
| `XeniumOvarianReproduction/scripts/plot_full_map_info_results.py` | **Current ovarian reporting** and the existing shared plotting-helper source imported by older dataset plotters; keep until helpers are moved to a neutral module. |
| `XeniumOvarianReproduction/scripts/plot_resolution_sweep.py` | **Current ovarian reporting**; keep. Its curve logic overlaps the E15.5 sweep plotters but its inputs/defaults are ovarian-specific. |
| `XeniumOvarianReproduction/scripts/plot_voronoi_clusters.py` | **Current ovarian reporting** and existing shared Voronoi-helper source; keep until helpers are moved to a neutral module. |

## Duplicated functionality and canonicalization decisions

### Label transfer

The clearest duplication is the E15.5/WholeEmbryo set:

```text
prepare_label_transfer.py
run_label_transfer_chunk.py
submit_label_transfer.py
check_label_transfer.py
combine_label_transfer.py
```

`run_label_transfer_chunk.py` differs only in its module docstring in the
audited copies. The WholeEmbryo versions of preparation, checking, and
combining contain important additional source-index, zero-variance, full-target
and transactional-output handling. Therefore the WholeEmbryo implementations
should be the behavioral basis for a future shared module, with expected sizes,
paths, section names, and UGE resources moved into task configuration. The old
`LabelTransfer/scripts/transfer_labels_knn.py` is not a duplicate of this
chunked workflow; it is a separate earlier one-shot STARmap pilot.

### Plotting

`LabelTransfer_E15.5/scripts/plot_mapinfo_sweep.py` and
`LabelTransfer_E15.5_WholeEmbryo/scripts/plot_mapinfo_sweep.py` share the same
curve/validation logic. The WholeEmbryo copy is the newer parameterized copy
and should become the existing-script canonical version until plotting helpers
are extracted.

The two E15.5 `plot_mapinfo_voronoi.py` files similarly share the same
Voronoi/score workflow. The WholeEmbryo copy is the newer superset because it
supports external labels, dataset titles, and source-row mappings; use it as
the behavioral basis for commonization.

The current ovarian `plot_full_map_info_results.py` and
`plot_voronoi_clusters.py` are already imported as helper modules by several
older dataset plotters. Their generic functions should eventually move to a
neutral shared plotting module; until then, do not delete the ovarian copies
just because older task folders are archived.

### Conversion

The conversion scripts are not interchangeable duplicates. STARmap RData,
GSE116470 DGE text, 10x Matrix Market, Qiu Matrix Market, and Visium HD HDF5
plus parquet each require different parsing and filtering. The repeated pieces
are validation idioms (`require`, sparse/count checks, gene-symbol collapsing,
AnnData round-trip verification), not whole programs. A future common helper
module may reduce this repetition, but no existing converter should be made the
single canonical converter for all datasets.

### Spatial preparation and job wrappers

The three families of spatial preparation scripts are intentionally different:
ovarian exports only Sections/XY; transferred-label preparation also exports
preexisting label arrays; ST-only and imputation preparation omit or alter
those labels. They should share validators eventually, but should remain
separate task entry points. The many `.sh` files are thin one-off resource/path
wrappers around either the core `max-info` CLI or dataset scripts; they are not
new algorithms and are the first candidates for consolidation.

## Original cleanup recommendations

1. Freeze provenance. Record which E15.5, WholeEmbryo, STARmap, and Visium HD
   results must remain reproducible, their source accessions, and the exact
   working-tree commit/checkpoint. Do not delete untracked task folders based
   only on their current Git status.
2. Add a small task manifest for each retained dataset containing raw inputs,
   processed outputs, config files, expected cell/gene counts, and the one
   supported command sequence. Use `LabelTransfer_E15.5/WORKFLOW.md` as the
   model.
3. Extract shared label-transfer code from the WholeEmbryo superset into a
   neutral module. Parameterize paths, expected dimensions, target mapping,
   resource requests, and section names; make both E15.5 entry points thin
   wrappers.
4. Extract shared plotting helpers from the ovarian scripts into a neutral
   plotting module. Keep dataset-specific CLI defaults in small wrappers and
   preserve the source-index/label-length checks.
5. Consolidate repeated AnnData validation helpers without merging the
   format-specific converters. Add syntax/unit tests for the shared helpers.
6. Replace one-off UGE shell wrappers with documented `max-info run` commands
   where the core CLI already provides the same behavior. Retain only wrappers
   that add a necessary resource request or a verified invariant check.
7. Archive, rather than immediately delete, the previous-pilot inspectors and
   plots, pilot jobs, one-off conversion jobs, and older E15.5 plotting copies.
   A proposed archive should be reviewed against result provenance first.
8. As a separate housekeeping change, remove generated `.DS_Store` and
   `__pycache__` artifacts and add/update ignore rules if needed. That cleanup
   was intentionally not done here.

The later reorganization implemented the archive-first portion of this plan
and extracted the current whole-embryo reference and shared label-transfer
scripts. See `MouseWholeEmbryo_MapInfo/README.md` for current paths and status.
