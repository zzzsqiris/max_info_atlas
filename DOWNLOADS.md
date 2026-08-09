# Removed downloaded datasets and generated artifacts

This inventory was written before cleanup on 2026-08-08. The active checkout is
the ovarian Xenium reproduction workflow; its source code, configs, raw input,
processed input, and results were kept. The removed material belonged to older
STARmap and E15.5 whole-embryo experiments and can be recreated from the
retained scripts and public accessions.

## Removed downloaded datasets

### GSE115746 ALM reference

- Removed files: `LabelTransfer/data/raw/GSE115746_cells_exon_counts.csv.gz`
  and the derived `LabelTransfer/data/processed/gse115746_alm.h5ad`.
- Original source: [NCBI GEO GSE115746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115746)
  (the file download URL was not recorded locally; the retained accession and
  `LabelTransfer/scripts/` conversion code identify the source file).
- Description: adult mouse anterior-lateral motor cortex single-cell exon-count
  reference used for the earlier STARmap label-transfer pilot.

### GSE116470 frontal-cortex reference

- Removed file: `LabelTransfer/data/raw/GSE116470_F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz`
  and the derived `LabelTransfer/data/processed/gse116470_frontal_cortex.h5ad`.
- Original download URL:
  <https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE116470&file=GSE116470_F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz&format=file>
- Description: processed adult mouse frontal-cortex Drop-seq matrix used as
  the scRNA reference for the earlier STARmap BZ9 pilot.

### E15.5 MERFISH spatial data

- Removed files under `LabelTransfer_E15.5/data/raw/`, including the two
  scRNA matrix bundles (`GSM6701497`, `GSM6701498`) and the MERFISH region-0
  cell-by-gene and cell-metadata files for `GSM7890130`.
- Original source: [NCBI GEO GSE247450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247450),
  sample `GSM7890130`.
- Description: E15.5 mouse whole-embryo MERFISH region used in the older
  label-transfer and Map Information comparisons.

### Qiu et al. whole-embryo scRNA reference

- Removed files under `LabelTransfer_E15.5_WholeEmbryo/data/raw/GSE228590/`,
  including the large gene-count matrix and per-run annotation tables.
- Original sources: [NCBI GEO GSE186069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186069)
  and [NCBI GEO GSE228590](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228590).
- Description: Qiu et al. whole-embryo single-cell reference and annotations
  used by the older E15.5 full-transfer, imputation, and spatial-only runs.

## Removed generated artifacts

- All generated data and results under `LabelTransfer/data/` and
  `LabelTransfer/results/`.
- All generated data and results under `LabelTransfer_E15.5/data/` and
  `LabelTransfer_E15.5/results/`.
- All generated data and results under
  `LabelTransfer_E15.5_WholeEmbryo/data/` and
  `LabelTransfer_E15.5_WholeEmbryo/results/`.
- Python bytecode caches and macOS `.DS_Store` files within those old
  experiment directories.

The corresponding source scripts, experiment configs, workflow documentation,
and the active `XeniumOvarianReproduction/` data/results were not removed.
