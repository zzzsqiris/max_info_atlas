# Mouse whole-embryo Map Information workflow

This directory is the unified experiment layer for whole-mouse-embryo Map
Information analyses. The reusable Map Information implementation remains in
`../src/max_info_atlases/`; this directory contains the shared reference,
spatial datasets, workflow wrappers, and historical experiments.

## Current design

```text
MouseWholeEmbryo_MapInfo/
├── reference/whole_embryo/       Nature 2024 E15.5 snRNA reference
├── datasets/visium_hd/           Paper 2; converted and ready for Hoffman2
├── datasets/stereo_seq/          Paper 1; downloaded, not yet inspected
├── shared/scripts/               Shared label-transfer and preparation code
└── archive/
    ├── starmap/                   Rejected: spatial dataset too small
    ├── merfish_ec_only/           Rejected: EC-only target has structural gaps
    └── merfish_whole_embryo/      Rejected: MERFISH panel is EC-biased
```

The active workflow for each retained spatial dataset is:

```text
whole-embryo snRNA reference
    -> reference clustering
    -> label transfer -> Map Information
    -> ST-only Map Information
    -> expression imputation -> Map Information
```

Dataset converters and thin wrappers remain dataset-specific. Shared
label-transfer workers are based on the later whole-embryo implementation,
including zero-variance filtering and source-index preservation. The files in
`archive/` are provenance snapshots and are not active pipeline entry points.

## Current status

- Visium HD conversion and local ST-only input preparation are complete.
- Visium HD Label Transfer is prepared for Hoffman2 but has not been executed.
- Stereo-seq source data is present locally under `data/raw/`, but has not yet
  been inspected, converted, or uploaded to Hoffman2.
- Visium HD imputation wrappers/configuration still need to be prepared.
- No historical data has been permanently deleted during this reorganization.

See `datasets/visium_hd/PAPER2_HOFFMAN_RUN.md` for the current Visium commands
and `HOFFMAN_REORGANIZE.md` for the one-time remote directory migration.
