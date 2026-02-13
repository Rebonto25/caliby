#!/usr/bin/env python3
"""
High-throughput sidechain rebuilding for many mutant PDBs.

Supports:
- FASPR (default, fastest) or SCWRL4
- Local parallel workers
- Optional SLURM array-style partitioning via stride
- Optional local scratch copy per file for faster shared-FS behavior

Example:
python rebuild_sidechains.py \
  --input-dir /path/to/mutant_pdbs_fast \
  --output-dir /path/to/mutant_pdbs_fullatom \
  --tool faspr \
  --workers 16 \
  --skip-existing
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Rebuild full-atom sidechains at scale.")
    p.add_argument("--input-dir", required=True, help="Input folder with PDB files")
    p.add_argument("--output-dir", required=True, help="Output folder for rebuilt PDB files")
    p.add_argument(
        "--tool",
        choices=["faspr", "scwrl4"],
        default="faspr",
        help="Rebuilder backend (faspr is fastest)",
    )
    p.add_argument("--tool-bin", default="", help="Optional explicit path to tool binary")
    p.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 1)),
        help="Parallel worker count",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip outputs that already exist and are non-empty",
    )
    p.add_argument(
        "--use-local-scratch",
        action="store_true",
        help="Copy each file to local scratch before running tool",
    )
    p.add_argument(
        "--scratch-dir",
        default="",
        help="Scratch base dir (default: $TMPDIR or /tmp)",
    )
    p.add_argument(
        "--manifest",
        default="",
        help="Optional text file (one PDB path per line). If omitted, scan input-dir.",
    )
    p.add_argument(
        "--array-id",
        type=int,
        default=None,
        help="Array task id. Default: $SLURM_ARRAY_TASK_ID or 0",
    )
    p.add_argument(
        "--array-count",
        type=int,
        default=None,
        help="Total number of array tasks. Default: $SLURM_ARRAY_TASK_COUNT or 1",
    )
    p.add_argument(
        "--failed-list",
        default="",
        help="Optional file to write failed input paths",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite outputs even if they exist",
    )
    return p.parse_args()


def _resolve_bin(tool: str, explicit: str) -> str:
    if explicit:
        bin_path = Path(explicit)
        if not bin_path.exists():
            raise FileNotFoundError(f"--tool-bin not found: {explicit}")
        return str(bin_path)

    candidates = []
    if tool == "faspr":
        candidates = ["FASPR", "faspr"]
    elif tool == "scwrl4":
        candidates = ["Scwrl4", "SCWRL4", "scwrl4"]

    for c in candidates:
        found = shutil.which(c)
        if found:
            return found
    raise FileNotFoundError(
        f"Could not find binary for tool='{tool}'. Add it to PATH or use --tool-bin."
    )


def _load_inputs(input_dir: Path, manifest: Optional[Path]) -> List[Path]:
    if manifest:
        files: List[Path] = []
        for line in manifest.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            p = Path(line)
            if p.is_file():
                files.append(p.resolve())
        return sorted(files)
    return sorted(p.resolve() for p in input_dir.glob("*.pdb") if p.is_file())


def _stride_partition(items: Sequence[Path], task_id: int, task_count: int) -> List[Path]:
    return [p for i, p in enumerate(items) if (i % task_count) == task_id]


def _tool_cmd(tool: str, tool_bin: str, in_pdb: Path, out_pdb: Path) -> List[str]:
    if tool == "faspr":
        return [tool_bin, "-i", str(in_pdb), "-o", str(out_pdb)]
    return [tool_bin, "-i", str(in_pdb), "-o", str(out_pdb)]


def _run_one(
    in_pdb: Path,
    output_dir: Path,
    tool: str,
    tool_bin: str,
    use_local_scratch: bool,
    scratch_dir: Path,
    skip_existing: bool,
    overwrite: bool,
) -> Tuple[Path, bool, str]:
    out_pdb = output_dir / in_pdb.name

    if out_pdb.exists() and out_pdb.stat().st_size > 0 and not overwrite:
        if skip_existing:
            return in_pdb, True, "skipped-existing"
        return in_pdb, False, "exists (use --overwrite or --skip-existing)"

    if use_local_scratch:
        scratch_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="rebuild_", dir=str(scratch_dir)) as td:
            td_path = Path(td)
            local_in = td_path / in_pdb.name
            local_out = td_path / f"out_{in_pdb.name}"
            shutil.copy2(in_pdb, local_in)
            cmd = _tool_cmd(tool, tool_bin, local_in, local_out)
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError:
                return in_pdb, False, "tool-failed"
            if not local_out.exists() or local_out.stat().st_size == 0:
                return in_pdb, False, "empty-output"
            shutil.move(str(local_out), str(out_pdb))
            return in_pdb, True, "ok"

    cmd = _tool_cmd(tool, tool_bin, in_pdb, out_pdb)
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return in_pdb, False, "tool-failed"

    if not out_pdb.exists() or out_pdb.stat().st_size == 0:
        return in_pdb, False, "empty-output"
    return in_pdb, True, "ok"


def _int_env(name: str, default: int) -> int:
    raw = os.environ.get(name, "")
    try:
        return int(raw)
    except ValueError:
        return default


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    manifest = Path(args.manifest).resolve() if args.manifest else None

    if not input_dir.is_dir() and manifest is None:
        raise FileNotFoundError(f"Input dir not found: {input_dir}")
    if manifest is not None and not manifest.is_file():
        raise FileNotFoundError(f"Manifest not found: {manifest}")

    output_dir.mkdir(parents=True, exist_ok=True)

    tool_bin = _resolve_bin(args.tool, args.tool_bin)

    array_id = args.array_id if args.array_id is not None else _int_env("SLURM_ARRAY_TASK_ID", 0)
    array_count = (
        args.array_count
        if args.array_count is not None
        else _int_env("SLURM_ARRAY_TASK_COUNT", 1)
    )
    if array_count < 1:
        raise ValueError("array-count must be >= 1")
    if not (0 <= array_id < array_count):
        raise ValueError(f"array-id must be in [0, {array_count - 1}]")

    all_inputs = _load_inputs(input_dir, manifest)
    if not all_inputs:
        print("No input PDB files found.")
        return

    task_inputs = _stride_partition(all_inputs, array_id, array_count)
    print(
        f"Found total={len(all_inputs)} PDBs; task {array_id}/{array_count} handles {len(task_inputs)}"
    )
    if not task_inputs:
        return

    workers = max(1, int(args.workers))
    workers = min(workers, len(task_inputs))

    scratch_base = Path(args.scratch_dir).resolve() if args.scratch_dir else Path(os.environ.get("TMPDIR", "/tmp")).resolve()

    results: List[Tuple[Path, bool, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [
            ex.submit(
                _run_one,
                p,
                output_dir,
                args.tool,
                tool_bin,
                bool(args.use_local_scratch),
                scratch_base,
                bool(args.skip_existing),
                bool(args.overwrite),
            )
            for p in task_inputs
        ]
        for fut in as_completed(futs):
            results.append(fut.result())

    ok = sum(1 for _, s, msg in results if s and msg != "skipped-existing")
    skipped = sum(1 for _, s, msg in results if s and msg == "skipped-existing")
    failed_items = [(p, msg) for p, s, msg in results if not s]

    if args.failed_list:
        failed_path = Path(args.failed_list).resolve()
    else:
        jid = os.environ.get("SLURM_ARRAY_JOB_ID", "nojid")
        tid = str(array_id)
        failed_path = output_dir / f"failed_{jid}_{tid}.txt"

    if failed_items:
        with failed_path.open("w") as f:
            for p, _ in failed_items:
                f.write(str(p) + "\n")

    print(f"Done. ok={ok}, skipped={skipped}, failed={len(failed_items)}")
    if failed_items:
        print(f"Failed list: {failed_path}")


if __name__ == "__main__":
    main()
