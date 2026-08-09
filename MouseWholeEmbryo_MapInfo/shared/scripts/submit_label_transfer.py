#!/usr/bin/env python3
"""Preview or submit a whole-embryo label-transfer UGE array job."""

import argparse
from pathlib import Path

from max_info_atlases.uge.submit import check_qsub_available, submit_array_job


REPO_DIR = Path(__file__).resolve().parents[3]
PIPELINE_DIR = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
VISIUM_RESULTS = PIPELINE_DIR / "datasets/visium_hd/results/label_transfer"
DEFAULT_CHUNKS_DIR = VISIUM_RESULTS / "work/chunks"
DEFAULT_WORK_DIR = VISIUM_RESULTS / "work"
DEFAULT_RESULTS_DIR = VISIUM_RESULTS / "chunk_results"
DEFAULT_LOGS_DIR = VISIUM_RESULTS / "logs"
DEFAULT_TEMP_DIR = VISIUM_RESULTS / "temp"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunks-dir", type=Path, default=DEFAULT_CHUNKS_DIR)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--logs-dir", type=Path, default=DEFAULT_LOGS_DIR)
    parser.add_argument("--temp-dir", type=Path, default=DEFAULT_TEMP_DIR)
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument("--job-name", default="e15_5_whole_label_transfer")
    parser.add_argument("--memory", default="8G")
    parser.add_argument("--runtime", default="12:00:00")
    parser.add_argument("--conda-env", default="xenium-mapinfo")
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Actually submit with qsub. Without this flag, print a dry run.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.k < 1:
        raise ValueError("--k must be positive.")
    if not args.chunks_dir.exists():
        raise FileNotFoundError(
            f"Missing chunks directory: {args.chunks_dir}. "
            "Run prepare_label_transfer.py first."
        )
    if args.submit and not check_qsub_available():
        raise RuntimeError("qsub is not available on this machine.")

    worker = (
        f'python "{SCRIPT_DIR / "run_label_transfer_chunk.py"}" '
        '--chunk-file "$chunk_file" '
        f'--work-dir "{args.work_dir}" '
        f'--output-dir "{args.results_dir}" '
        f'--k {args.k}'
    )
    job_id, _ = submit_array_job(
        chunks_dir=args.chunks_dir,
        base_dir=REPO_DIR,
        temp_dir=args.temp_dir,
        logs_dir=args.logs_dir,
        job_name=args.job_name,
        memory=args.memory,
        runtime=args.runtime,
        conda_env=args.conda_env,
        worker_command=worker,
        dry_run=not args.submit,
    )
    if job_id:
        print(f"label_transfer_job_id: {job_id}")


if __name__ == "__main__":
    main()
