# max_info_atlas project overview

Last updated: 2026-08-08

This is the persistent orientation note for future work in
`/Users/shiqizhong/Code/max_info_atlas`. It describes the current local layout,
the matching Hoffman2 layout, the canonical Map Information implementation,
and the status of the whole-mouse-embryo workflows.

## Repository and Git status

- Local repository: `/Users/shiqizhong/Code/max_info_atlas`
- Hoffman2 repository: `~/project-rwollman/max_info_atlas`
- User GitHub remote (`origin`): `https://github.com/zzzsqiris/max_info_atlas.git`
- Professor repository (`upstream`): `https://github.com/rwollman/max_info_atlas.git`
- Current local branch at this update: `ovarian-xenium-reproduction`
- Recommended branch for this reorganization: `mouse-whole-embryo-mapinfo`

The whole-embryo reorganization is present in the working tree but has not yet
been claimed as committed or pushed. Use GitHub Desktop to create the new
branch and commit only the whole-embryo files. The following pre-existing
working-tree changes are unrelated and must not be included in that commit:

```text
XeniumOvarianReproduction/scripts/plot_full_map_info_results.py
XeniumOvarianReproduction/scripts/plot_voronoi_clusters.py
src/max_info_atlases/features/celltype.py
```

The `.gitignore` excludes all whole-embryo `data/`, `results/`, `logs/`,
`figures/`, and `tmp/` directories. The full historical experiment snapshots
under `MouseWholeEmbryo_MapInfo/archive/` also remain local/Hoffman2-only;
GitHub tracks only `archive/README.md`. Raw or processed `.h5ad`, `.npy`, HDF5,
parquet, compressed archives, and generated outputs must not be committed.

## Current repository structure

```text
max_info_atlas/
├── src/max_info_atlases/             Canonical reusable package and CLI
├── config/                           Core and ovarian run configurations
├── tests/                            Unit and small pipeline tests
├── XeniumOvarianReproduction/       Separate ovarian Xenium workflow
├── MouseWholeEmbryo_MapInfo/        Unified whole-embryo experiment layer
│   ├── reference/
│   │   └── whole_embryo/             Nature 2024 E15.5 snRNA reference
│   ├── datasets/
│   │   ├── visium_hd/                Paper 2, current active spatial target
│   │   └── stereo_seq/               Paper 1, downloaded but not inspected
│   ├── shared/
│   │   └── scripts/                  Shared Label Transfer/ST preparation
│   └── archive/
│       ├── starmap/                  Superseded pilot
│       ├── merfish_ec_only/          Superseded EC-only workflow
│       └── merfish_whole_embryo/     Superseded MERFISH target workflow
├── DOWNLOADS.md                      Historical download/provenance notes
├── README.md                         Package installation and CLI overview
├── TESTING.md                        Core testing instructions
└── PROJECT_OVERVIEW.md               This document
```

The dataset folders are experiment, configuration, and provenance layers.
They do not replace the reusable implementation under `src/max_info_atlases/`.
Prefer dataset-specific wrappers/configuration over modifying `src/`.

## Canonical Map Information pipeline

The reusable package is exposed through the `max-info` command configured in
`setup.py`:

```text
AnnData + Sections.npy + <section>_XY.npy
    -> features
    -> graphs
    -> clustering
    -> percolation
    -> aggregation / reduced score CSV
```

Normal configured commands are:

```bash
max-info run summary --config <config.yaml>
max-info run prepare --config <config.yaml>
max-info run submit --config <config.yaml> --step <step>
max-info run status --config <config.yaml>
```

Current configurations generally use PCA50, correlation kNN with `k=15`, a
50-point Leiden resolution sweep, `max_k=500`, `pbond_steps=101`, and Leiden
`random_seed: 42`. Do not treat older runs without the explicit seed as
identical deterministic results.

## Unified whole-embryo workflow

The intended workflow for each retained spatial dataset is:

```text
Nature 2024 whole-embryo snRNA reference
    -> reference clustering
    -> Label Transfer -> Map Information
    -> ST-only Map Information
    -> expression imputation -> Map Information
```

The later WholeEmbryo implementation is the behavioral basis for shared Label
Transfer code because it preserves source-row mappings, excludes zero-variance
shared-expression rows safely, represents excluded target cells explicitly,
checks chunk completeness, and writes combined outputs transactionally.

Primary transferred labels use majority voting across the `k=15` reference
neighbors. Weighted labels, vote confidence, and neighbor-distance summaries
are diagnostics unless a dataset-specific analysis explicitly chooses them.

### Shared active scripts

```text
MouseWholeEmbryo_MapInfo/shared/scripts/
├── prepare_label_transfer.py
├── run_label_transfer_chunk.py
├── submit_label_transfer.py
├── check_label_transfer.py
├── combine_label_transfer.py
└── prepare_st_only_full_mapinfo.py
```

These active copies have been updated for the unified directory layout. Older
STARmap/MERFISH-specific versions remain only in the local/Hoffman2 archive.
Imputation scripts were not promoted into the active shared area because their
existing defaults and invariants remain MERFISH-specific; they must be adapted
separately for Visium HD and Stereo-seq.

## Shared whole-embryo reference

Active paths:

```text
MouseWholeEmbryo_MapInfo/reference/whole_embryo/data/processed/
├── qiu_e15_5_whole_embryo_reference_raw.h5ad
└── qiu_e15_5_whole_embryo_reference_sections.npy
```

Verified properties:

- Reference cells: 215,284
- Section: `Qiu_E15.5_embryo_59`
- Local reference directory size: approximately 3.3 GB
- The AnnData was successfully opened read-only before reorganization.
- File sizes and HDF5 signature were rechecked after the local move.

Hoffman2 reference clustering was moved to:

```text
MouseWholeEmbryo_MapInfo/reference/whole_embryo/results/reference_clustering/
```

The user-run Hoffman2 verification found exactly 50
`Qiu_E15.5_embryo_59.npy` reference label arrays. Reuse this completed
reference clustering; do not rerun it unless those outputs later become
missing or invalid.

## Paper 2: Visium HD

Active directory:

```text
MouseWholeEmbryo_MapInfo/datasets/visium_hd/
```

Current validated inputs:

- Converted AnnData: `data/processed/mouse_embryo_visium_hd_raw.h5ad`
- Shape: 449,700 spatial observations × 33,670 genes
- Sparse CSR expression matrix
- ST-only sections: `data/processed/visium_hd_st_only_mapinfo/Sections.npy`
- Coordinates:
  `data/processed/visium_hd_st_only_mapinfo/xy_coordinates/MouseEmbryo_VisiumHD_XY.npy`
- Coordinate shape: 449,700 × 2
- Section: `MouseEmbryo_VisiumHD`
- Total local Visium directory size, including raw archives: approximately
  7.0 GB

Current readiness:

- Raw-data conversion: complete and validated
- Local ST-only spatial-input preparation: complete
- Hoffman2 data upload: complete
- Hoffman2 code/config synchronization: complete
- Hoffman2 active Python/shell syntax check: passed under Python 3.11.15
- Reference clustering: complete and reusable (50 label arrays)
- Label Transfer: ready to start on Hoffman2, but not yet executed
- ST-only Map Information: configuration and inputs ready, but not executed
- Imputation + Map Information: not yet adapted or ready

The complete current run instructions are in:

```text
MouseWholeEmbryo_MapInfo/datasets/visium_hd/PAPER2_HOFFMAN_RUN.md
```

The next Label Transfer action, when explicitly requested, is the Visium
preparation qsub job. Do not rerun conversion or reference clustering first.

## Paper 1: Stereo-seq / MOSTA

Active directory:

```text
MouseWholeEmbryo_MapInfo/datasets/stereo_seq/
```

Current local source:

```text
data/raw/Mouse_embryo_all_stage.h5ad
```

Current status:

- Download is present locally and occupies approximately 21 GB.
- The file was moved without renaming into the standard `data/raw/` location.
- Its contents have not yet been inspected in this workflow.
- No conversion, Label Transfer preparation, ST-only preparation, or
  imputation adaptation has been performed.
- The source has not been uploaded to Hoffman2.

Do not infer its AnnData schema, embryo stages, spatial coordinate fields, or
pipeline readiness until a separate lightweight inspection is requested.

## Historical workflows

The rejected workflows are retained for provenance but are not active entry
points:

| Archive | Reason superseded | Role of retained material |
|---|---|---|
| `archive/starmap/` | Spatial dataset was too small | Pilot converters, one-shot transfer code, configs, and plotting provenance |
| `archive/merfish_ec_only/` | EC-only target produced structural gaps | Earlier chunked Label Transfer, ST-only, imputation, configs, and results provenance |
| `archive/merfish_whole_embryo/` | MERFISH gene panel was biased toward endothelial cells | Later WholeEmbryo implementation used as the behavioral source for shared code |

The full archive snapshots remain local and on Hoffman2 but are excluded from
GitHub. No historical analysis data was permanently deleted during the
reorganization. Generated `.DS_Store` and `__pycache__` files and redundant
active copies of unadapted MERFISH scripts were removed; canonical originals
remain in the archive.

## Local and Hoffman2 synchronization

Both locations now use the unified `MouseWholeEmbryo_MapInfo/` hierarchy.
The one-time migration and sync record is:

```text
MouseWholeEmbryo_MapInfo/HOFFMAN_REORGANIZE.md
```

The user-run Hoffman2 checks verified:

- all five required reference/Visium input files exist;
- file sizes are consistent with the local copies;
- 50 reference clustering label arrays exist;
- synchronized active Python scripts parse successfully in the
  `xenium-mapinfo` environment;
- active shell wrappers pass `bash -n`.

Hoffman2 jobs should run from:

```bash
cd ~/project-rwollman/max_info_atlas
. /u/local/Modules/default/init/modules.sh
module load conda
conda activate xenium-mapinfo
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
```

The login-node default `python` is Python 2 and does not provide `pathlib`.
Always activate `xenium-mapinfo` before Python validation or workflow commands.

## Data and Git rules

- Never commit raw/processed biological data or generated pipeline outputs.
- Track scripts, configs, README/workflow documents, and small provenance
  manifests only.
- Do not use GitHub Desktop's default "select all" if unrelated working-tree
  edits are present.
- Push the whole-embryo work to the user's `origin`, not the professor's
  `upstream`.
- Do not claim a workflow is committed or published until Git status and the
  remote branch confirm it.
- Do not claim a submitted HPC job succeeded from its job ID alone; validate
  outputs, manifests, logs, score counts, and row/coordinate alignment.
- Keep AnnData row order as the primary alignment key and preserve source-index
  arrays whenever zero-variance cells are removed from Label Transfer.

## Current next steps

1. In GitHub Desktop, create `mouse-whole-embryo-mapinfo`, commit only the
   whole-embryo reorganization plus `.gitignore`, `DOWNLOADS.md`, and this
   overview, then publish to `origin`.
2. When execution is requested, start Paper 2 Visium HD Label Transfer on
   Hoffman2 using `PAPER2_HOFFMAN_RUN.md`; reuse the existing reference labels.
3. Prepare and run Paper 2 ST-only Map Information separately.
4. Inspect Paper 1 Stereo-seq before designing its converter or wrappers.
5. Adapt imputation for each retained spatial platform only after the earlier
   stages are ready; do not reuse MERFISH-specific fixed counts silently.
