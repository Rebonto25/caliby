#!/usr/bin/env python3
"""
Simple sidechain rebuilding with FASPR.

Usage:
  # Single file
  python rebuild_sidechains_simple.py -i input.pdb -o output.pdb

  # Directory (single-threaded)
  python rebuild_sidechains_simple.py --input-dir ./pdbs --output-dir ./out

  # Directory (parallel)
  python rebuild_sidechains_simple.py --input-dir ./pdbs --output-dir ./out --workers 8

  # With explicit FASPR path
  python rebuild_sidechains_simple.py --input-dir ./pdbs --output-dir ./out --faspr /path/to/FASPR
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def find_faspr(explicit: str | None = None) -> Path:
    """Locate FASPR binary."""
    if explicit:
        p = Path(explicit)
        if p.is_file():
            return p
        raise FileNotFoundError(f"FASPR not found at: {explicit}")

    for name in ("FASPR", "faspr"):
        found = shutil.which(name)
        if found:
            return Path(found)

    raise FileNotFoundError(
        "FASPR not found in PATH. Use --faspr to specify location."
    )


def run_faspr(faspr_bin: Path, in_pdb: Path, out_pdb: Path, seq: str | None = None) -> tuple[Path, bool, str]:
    """
    Run FASPR on a single PDB.

    Returns (input_path, success, message).
    """
    cmd = [str(faspr_bin), "-i", str(in_pdb), "-o", str(out_pdb)]

    # Optional: provide sequence for mutations or fixing residues
    if seq:
        seq_file = out_pdb.with_suffix(".seq.txt")
        seq_file.write_text(seq)
        cmd.extend(["-s", str(seq_file)])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            return in_pdb, False, f"FASPR error: {result.stderr.strip()}"
        if not out_pdb.exists() or out_pdb.stat().st_size == 0:
            return in_pdb, False, "empty output"
        return in_pdb, True, "ok"
    except subprocess.TimeoutExpired:
        return in_pdb, False, "timeout"
    except Exception as e:
        return in_pdb, False, str(e)


def process_single(args: argparse.Namespace) -> int:
    """Process a single input/output pair."""
    faspr = find_faspr(args.faspr)
    in_pdb = Path(args.input).resolve()
    out_pdb = Path(args.output).resolve()

    if not in_pdb.exists():
        print(f"ERROR: input not found: {in_pdb}")
        return 1

    out_pdb.parent.mkdir(parents=True, exist_ok=True)

    _, ok, msg = run_faspr(faspr, in_pdb, out_pdb, args.seq)
    if ok:
        print(f"OK: {out_pdb}")
        return 0
    else:
        print(f"FAILED: {in_pdb} ({msg})")
        return 1


def process_directory(args: argparse.Namespace) -> int:
    """Process all PDBs in a directory."""
    faspr = find_faspr(args.faspr)
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    if not input_dir.is_dir():
        print(f"ERROR: input directory not found: {input_dir}")
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    # Gather inputs
    pdbs = sorted(input_dir.glob("*.pdb"))
    if not pdbs:
        print(f"No PDB files found in {input_dir}")
        return 0

    # Filter existing if requested
    if args.skip_existing:
        pdbs = [p for p in pdbs if not (output_dir / p.name).exists()]

    if not pdbs:
        print("All outputs already exist (--skip-existing).")
        return 0

    print(f"Processing {len(pdbs)} PDB files...")

    workers = max(1, args.workers)
    ok_count = 0
    fail_count = 0
    failed_paths: list[Path] = []

    if workers == 1:
        # Single-threaded
        for pdb in pdbs:
            out_pdb = output_dir / pdb.name
            _, ok, msg = run_faspr(faspr, pdb, out_pdb)
            if ok:
                ok_count += 1
            else:
                fail_count += 1
                failed_paths.append(pdb)
                if args.verbose:
                    print(f"FAILED: {pdb.name} ({msg})")
    else:
        # Parallel
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(run_faspr, faspr, pdb, output_dir / pdb.name): pdb
                for pdb in pdbs
            }
            for fut in as_completed(futures):
                in_pdb, ok, msg = fut.result()
                if ok:
                    ok_count += 1
                else:
                    fail_count += 1
                    failed_paths.append(in_pdb)
                    if args.verbose:
                        print(f"FAILED: {in_pdb.name} ({msg})")

    print(f"Done. ok={ok_count}, failed={fail_count}")

    if failed_paths and args.failed_list:
        failed_file = Path(args.failed_list)
        failed_file.write_text("\n".join(str(p) for p in failed_paths) + "\n")
        print(f"Failed list: {failed_file}")

    return 0 if fail_count == 0 else 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rebuild sidechains with FASPR.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Single file mode
    parser.add_argument("-i", "--input", help="Single input PDB file")
    parser.add_argument("-o", "--output", help="Single output PDB file")
    parser.add_argument("-s", "--seq", help="Optional sequence string for FASPR -s flag")

    # Directory mode
    parser.add_argument("--input-dir", help="Input directory with PDB files")
    parser.add_argument("--output-dir", help="Output directory for rebuilt PDBs")

    # Common options
    parser.add_argument("--faspr", help="Path to FASPR binary (auto-detected if in PATH)")
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers (default: 1)")
    parser.add_argument("--skip-existing", action="store_true", help="Skip existing outputs")
    parser.add_argument("--verbose", action="store_true", help="Print failures as they happen")
    parser.add_argument("--failed-list", help="File to write failed input paths")

    args = parser.parse_args()

    # Decide mode
    if args.input and args.output:
        return process_single(args)
    elif args.input_dir and args.output_dir:
        return process_directory(args)
    else:
        parser.print_help()
        print("\nERROR: Provide either (-i, -o) or (--input-dir, --output-dir)")
        return 1


if __name__ == "__main__":
    sys.exit(main())
