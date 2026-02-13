#!/usr/bin/env python3
"""
Generate mutated PDBs from a parent antibody-antigen complex PDB and FASTA sequences.

For each FASTA entry, this script mutates a single chain (default: H) by residue-name
substitution only (coordinates are kept from the parent structure). The rest of the
complex (antigen + other chains) is preserved as-is.

Key features:
- Checks sequence validity and chain/length compatibility
- Checks mutation count (default max: 17)
- Optional multiprocessing
- Reuses loaded parent structure per process (avoids re-reading PDB for every output)
"""

from __future__ import annotations

import argparse
import copy
import csv
import os
import re
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

# Global process-local cache (initialized once per process).
_WORKER_STATE: Dict[str, object] = {}


THREE_TO_ONE = {
	"ALA": "A",
	"ARG": "R",
	"ASN": "N",
	"ASP": "D",
	"CYS": "C",
	"GLN": "Q",
	"GLU": "E",
	"GLY": "G",
	"HIS": "H",
	"ILE": "I",
	"LEU": "L",
	"LYS": "K",
	"MET": "M",
	"PHE": "F",
	"PRO": "P",
	"SER": "S",
	"THR": "T",
	"TRP": "W",
	"TYR": "Y",
	"VAL": "V",
	"MSE": "M",  # common selenium methionine
}

ONE_TO_THREE = {
	"A": "ALA",
	"R": "ARG",
	"N": "ASN",
	"D": "ASP",
	"C": "CYS",
	"Q": "GLN",
	"E": "GLU",
	"G": "GLY",
	"H": "HIS",
	"I": "ILE",
	"L": "LEU",
	"K": "LYS",
	"M": "MET",
	"F": "PHE",
	"P": "PRO",
	"S": "SER",
	"T": "THR",
	"W": "TRP",
	"Y": "TYR",
	"V": "VAL",
}

VALID_AA = set(ONE_TO_THREE.keys())
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


@dataclass(frozen=True)
class Job:
	seq_id: str
	sequence: str
	out_path: str
	max_mutations: int
	keep_sidechain: bool


def _normalize_sequence(seq: str) -> str:
	seq = seq.upper().replace(" ", "").replace("\n", "")
	seq = seq.replace("-", "").replace(".", "")
	return seq


def _safe_stem(text: str) -> str:
	stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
	return stem or "variant"


def _require_biopython():
	try:
		from Bio import SeqIO  # type: ignore
		from Bio.PDB import PDBIO, PDBParser  # type: ignore
		from Bio.PDB.Polypeptide import is_aa  # type: ignore
	except ImportError as e:
		raise ModuleNotFoundError(
			"Biopython is required. Install with: pip install biopython"
		) from e
	return SeqIO, PDBIO, PDBParser, is_aa


def _protein_residues(chain) -> List:
	is_aa_func = _WORKER_STATE.get("is_aa_func")
	if is_aa_func is None:
		_, _, _, is_aa_func = _require_biopython()
	residues = []
	for res in chain:
		if res.id[0] == " " and is_aa_func(res, standard=False):
			residues.append(res)
	return residues


def _chain_sequence(chain) -> str:
	seq = []
	for res in _protein_residues(chain):
		aa = THREE_TO_ONE.get(res.resname.upper())
		if aa is None:
			aa = "X"
		seq.append(aa)
	return "".join(seq)


def _mutation_list(parent_seq: str, mutant_seq: str) -> List[Tuple[int, str, str]]:
	return [
		(i + 1, p, m)
		for i, (p, m) in enumerate(zip(parent_seq, mutant_seq))
		if p != m
	]


def _trim_mutated_sidechains(residue, new_aa: str) -> None:
	"""Keep only backbone atoms for mutated residues (+CB for non-GLY)."""
	keep = set(BACKBONE_ATOMS)
	if new_aa != "G":
		keep.add("CB")
	for atom in list(residue.get_atoms()):
		if atom.id not in keep:
			residue.detach_child(atom.id)


def _init_worker(parent_pdb: str, chain_id: str) -> None:
	_, PDBIO, PDBParser, is_aa_func = _require_biopython()
	parser = PDBParser(QUIET=True)
	structure = parser.get_structure("parent", parent_pdb)
	model = structure[0]
	if chain_id not in model:
		raise ValueError(f"Chain '{chain_id}' not found in parent PDB: {parent_pdb}")

	chain = model[chain_id]
	seq = _chain_sequence(chain)
	if not seq:
		raise ValueError(f"No protein residues found in chain '{chain_id}'.")
	if "X" in seq:
		raise ValueError(
			f"Chain '{chain_id}' has unsupported residues in parent PDB (mapped to X)."
		)

	_WORKER_STATE["template_structure"] = structure
	_WORKER_STATE["chain_id"] = chain_id
	_WORKER_STATE["parent_seq"] = seq
	_WORKER_STATE["PDBIO"] = PDBIO
	_WORKER_STATE["is_aa_func"] = is_aa_func


def _run_job(job: Job) -> Dict[str, object]:
	parent_seq = _WORKER_STATE["parent_seq"]
	chain_id = _WORKER_STATE["chain_id"]
	template = _WORKER_STATE["template_structure"]

	mutant_seq = _normalize_sequence(job.sequence)

	invalid = sorted({c for c in mutant_seq if c not in VALID_AA})
	if invalid:
		return {
			"id": job.seq_id,
			"status": "failed",
			"reason": f"invalid AA letters: {''.join(invalid)}",
			"mutations": "",
			"num_mutations": 0,
			"out_path": job.out_path,
		}

	if len(mutant_seq) != len(parent_seq):
		return {
			"id": job.seq_id,
			"status": "failed",
			"reason": (
				f"length mismatch (parent={len(parent_seq)}, mutant={len(mutant_seq)})"
			),
			"mutations": "",
			"num_mutations": 0,
			"out_path": job.out_path,
		}

	muts = _mutation_list(parent_seq, mutant_seq)
	if len(muts) > job.max_mutations:
		return {
			"id": job.seq_id,
			"status": "failed",
			"reason": (
				f"too many mutations: {len(muts)} > max {job.max_mutations}"
			),
			"mutations": ";".join(f"{p}{i}{m}" for i, p, m in muts),
			"num_mutations": len(muts),
			"out_path": job.out_path,
		}

	structure = copy.deepcopy(template)
	chain = structure[0][chain_id]
	residues = _protein_residues(chain)

	for idx_1based, _old, new in muts:
		res = residues[idx_1based - 1]
		res.resname = ONE_TO_THREE[new]
		if not job.keep_sidechain:
			_trim_mutated_sidechains(res, new)

	pdbio_cls = _WORKER_STATE.get("PDBIO")
	if pdbio_cls is None:
		_, pdbio_cls, _, _ = _require_biopython()
	io = pdbio_cls()
	io.set_structure(structure)
	io.save(job.out_path)

	mut_str = ";".join(f"{p}{i}{m}" for i, p, m in muts)
	reason = "ok"
	if len(muts) == 0:
		reason = "no differences vs parent"

	return {
		"id": job.seq_id,
		"status": "ok",
		"reason": reason,
		"mutations": mut_str,
		"num_mutations": len(muts),
		"out_path": job.out_path,
	}


def _iter_fasta_jobs(
	fasta_path: Path,
	out_dir: Path,
	max_mutations: int,
	keep_sidechain: bool,
	overwrite: bool,
) -> List[Job]:
	SeqIO, _, _, _ = _require_biopython()
	jobs: List[Job] = []
	seen: Dict[str, int] = {}

	for rec in SeqIO.parse(str(fasta_path), "fasta"):
		base = _safe_stem(rec.id)
		seen[base] = seen.get(base, 0) + 1
		seq_id = base if seen[base] == 1 else f"{base}_{seen[base]}"

		out_path = out_dir / f"{seq_id}.pdb"
		if out_path.exists() and not overwrite:
			continue

		jobs.append(
			Job(
				seq_id=seq_id,
				sequence=str(rec.seq),
				out_path=str(out_path),
				max_mutations=max_mutations,
				keep_sidechain=keep_sidechain,
			)
		)

	return jobs


def parse_args() -> argparse.Namespace:
	p = argparse.ArgumentParser(
		description="Generate mutated PDBs by substituting one chain sequence from FASTA."
	)
	p.add_argument("--parent-pdb", required=True, help="Parent complex PDB file")
	p.add_argument("--fasta", required=True, help="FASTA of mutant chain sequences")
	p.add_argument("--out-dir", required=True, help="Output folder for mutated PDBs")
	p.add_argument("--chain", default="H", help="Chain ID to mutate (default: H)")
	p.add_argument(
		"--max-mutations",
		type=int,
		default=17,
		help="Maximum allowed mutations vs parent chain sequence (default: 17)",
	)
	p.add_argument(
		"--processes",
		type=int,
		default=os.cpu_count() or 1,
		help="Number of processes (default: all available cores)",
	)
	p.add_argument(
		"--keep-sidechain",
		action="store_true",
		help=(
			"Keep all original sidechain atoms in mutated residues. "
			"By default only backbone(+CB for non-GLY) is kept for mutated positions."
		),
	)
	p.add_argument(
		"--overwrite",
		action="store_true",
		help="Overwrite existing output PDB files",
	)
	p.add_argument(
		"--summary-csv",
		default="summary.csv",
		help="Summary CSV filename written into --out-dir (default: summary.csv)",
	)
	return p.parse_args()


def _write_summary(path: Path, rows: Iterable[Dict[str, object]]) -> None:
	fields = ["id", "status", "reason", "num_mutations", "mutations", "out_path"]
	with path.open("w", newline="") as f:
		w = csv.DictWriter(f, fieldnames=fields)
		w.writeheader()
		for r in rows:
			w.writerow({k: r.get(k, "") for k in fields})


def main() -> None:
	args = parse_args()

	parent_pdb = Path(args.parent_pdb).resolve()
	fasta_path = Path(args.fasta).resolve()
	out_dir = Path(args.out_dir).resolve()

	if not parent_pdb.is_file():
		raise FileNotFoundError(f"Parent PDB not found: {parent_pdb}")
	if not fasta_path.is_file():
		raise FileNotFoundError(f"FASTA not found: {fasta_path}")
	if len(args.chain) != 1:
		raise ValueError("--chain must be a single character chain ID")
	if args.max_mutations < 0:
		raise ValueError("--max-mutations must be >= 0")

	out_dir.mkdir(parents=True, exist_ok=True)

	jobs = _iter_fasta_jobs(
		fasta_path=fasta_path,
		out_dir=out_dir,
		max_mutations=args.max_mutations,
		keep_sidechain=args.keep_sidechain,
		overwrite=args.overwrite,
	)
	if not jobs:
		print("No jobs to run (all outputs may already exist).")
		return

	n_proc = max(1, int(args.processes))
	n_proc = min(n_proc, len(jobs))

	results: List[Dict[str, object]] = []
	if n_proc == 1:
		_init_worker(str(parent_pdb), args.chain)
		for j in jobs:
			results.append(_run_job(j))
	else:
		with ProcessPoolExecutor(
			max_workers=n_proc,
			initializer=_init_worker,
			initargs=(str(parent_pdb), args.chain),
		) as ex:
			for out in ex.map(_run_job, jobs):
				results.append(out)

	summary_path = out_dir / args.summary_csv
	_write_summary(summary_path, results)

	ok = sum(1 for r in results if r["status"] == "ok")
	failed = len(results) - ok
	print(f"Done. ok={ok}, failed={failed}")
	print(f"Summary: {summary_path}")

	if failed:
		print("Failed entries:")
		for r in results:
			if r["status"] != "ok":
				print(f"  - {r['id']}: {r['reason']}")


if __name__ == "__main__":
	main()
