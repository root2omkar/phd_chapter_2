#!/usr/bin/env python3

from __future__ import annotations
from typing import List, Set, Tuple, Dict
from pathlib import Path
from Bio import SeqIO
from itertools import combinations
import gzip, argparse, sys, re, math

# ---------------------------
# Auxiliary functions
# ---------------------------

# Parse CLI arguments: FASTA files (-f/--fastas), verbosity (-v), and reporting options.

def get_args():
	"""Reads different alignments in FASTA format and compares them."""
	parser = argparse.ArgumentParser(
		description="Reads different alignments in FASTA format and compares them."
	)
	parser.add_argument(
		"-i", "--input",
		type=str,
		nargs='+',
		help="List FASTA files separated by space",
		required = True
	)
	parser.add_argument(
		"-v", "--verbose",
		action="store_true",
		help="Increase output verbosity"
	)
	return parser.parse_args()

# Print to stdout (normal progress/info).

def po(message: str) -> None:
	sys.stdout.write(f"{message}\n")
	return

# Print to stderr (warnings/errors/verbose diagnostics).

def pe(message: str) -> None:
	sys.stderr.write(f"{message}\n")
	return

# Normalize a raw sequence: uppercase, remove spaces, standardize ambiguous symbols if needed.

def cleanSeq(seq: str) -> str:
	clean_seq = re.sub(r"[-–—]+", "", seq)
	clean_seq = re.sub(r"^[Nn\?]+", "", clean_seq)
	clean_seq = re.sub(r"[Nn\?]+$", "", clean_seq)
	return clean_seq

# Convert terminal N's to '-' (gap-like) or otherwise harmonize dash handling before comparisons.

def handleDashes(seq: str) -> str:
	new_seq = re.sub(r"[-–—]", "-", seq)
	new_seq = re.sub(r'^N+', lambda m: '-' * len(m.group(0)), new_seq)
	new_seq = re.sub(r'N+$', lambda m: '-' * len(m.group(0)), new_seq)
	return new_seq

# Encode a sequence into per-position keys like 'seqID:colIndex' used to build position-wise sets.

def seq2pos(seqid: str, seq: str) -> List[str]:
	list_positions: List[str] = []
	position = 0
	for char in seq:
		if char == "-":
			code = f"{seqid}:-"
		else:
			code = f"{seqid}:{position}"
			position += 1
		list_positions.append(code)
	return list_positions

# Encode a sequence into per-character-state labels like 'seqID:STATE' used to build state-wise sets.

def seq2char(seqid: str, seq: str) -> List[str]:
	list_states: List[str] = []
	for char in seq:
		code = f"{seqid}:{char}"
		list_states.append(code)
	return list_states

# Jaccard similarity between two sets: |A∩B| / |A∪B|; robust for comparing set-of-sets.

def jaccard(a: Set[str], b: Set[str]) -> float:
	if not a and not b:
		return 1.0
	inter = len(a & b)
	union = len(a | b)
	return inter / union if union else 0.0

def hungarian(cost):
		n = len(cost)
		u = [0.0]*(n+1)
		v = [0.0]*(n+1)
		p = [0]*(n+1)
		way = [0]*(n+1)
		for i in range(1, n+1):
			p[0] = i
			j0 = 0
			minv = [math.inf]*(n+1)
			used = [False]*(n+1)
			while True:
				used[j0] = True
				i0 = p[j0]
				delta = math.inf
				j1 = 0
				for j in range(1, n+1):
					if not used[j]:
						cur = cost[i0-1][j-1] - u[i0] - v[j]
						if cur < minv[j]:
							minv[j] = cur; way[j] = j0
						if minv[j] < delta:
							delta = minv[j]; j1 = j
				for j in range(n+1):
					if used[j]:
						u[p[j]] += delta; v[j] -= delta
					else:
						minv[j] -= delta
				j0 = j1
				if p[j0] == 0:
					break
			while True:
				j1 = way[j0]
				p[j0] = p[j1]
				j0 = j1
				if j0 == 0:
					break
		match = [0]*(n+1)
		for j in range(1, n+1):
			match[p[j]] = j
		# match[i] = assigned column for row i (1-indexed to 1-indexed)
		return match[1:]  # drop dummy 0

# Solve a one-to-one matching maximizing total Jaccard similarity (Hungarian algorithm under the hood).

def optimal_match_score(U: List[Set[str]], V: List[Set[str]]) -> Tuple[float, list]:
	"""
	Returns (score, matching) where score in [0,1] and matching is list of (i,j,similarity).
	Uses a simple Hungarian implementation (square matrix required).
	Example:
	U = [ {"a","b"}, {"c"}, {"x","y"} ]
	V = [ {"b","a"}, {"z"}, {"x"} ]
	score, matching = optimal_match_score(U, V)
	print(score, matching)
	"""
	# Build similarity matrix
	m, n = len(U), len(V)
	k = max(m, n)
	S = [[0.0]*k for _ in range(k)]
	for i in range(m):
		for j in range(n):
			S[i][j] = jaccard(U[i], V[j])
	# Convert to cost = 1 - similarity
	C = [[1.0 - S[i][j] for j in range(k)] for i in range(k)]
	# Hungarian algorithm (O(k^3)), pure Python
	# Adapted from a standard minimal implementation for teaching purposes
	# — fine for class sizes up to a few hundred.
	assignment = hungarian(C)  # list of length k; assignment[i-1] in 1..k
	# Build matching and compute normalized score
	pairs = []
	sim_sum = 0.0
	for i in range(k):
		j = assignment[i] - 1
		if i < m and j < n:
			sim = S[i][j]
			pairs.append((i, j, sim))
			sim_sum += sim
		else:
			# matched to a dummy => similarity 0
			pairs.append((i, None if j >= n else j, 0.0))
	score = sim_sum / k if k else 1.0
	return score, pairs

# Combine exact matches and optimal matching of unique sets into a single similarity score in [0,1].

def final_similarity(shared_count: int, m: int, n: int, match_score: float) -> float:
	k = max(m, n)
	total = shared_count + m + n
	if total == 0:
		return 1.0  # both lists empty → perfect similarity
	return (float(shared_count) + match_score * float(k)) / float(total)

# ---------------------------
# Class Alignments
# ---------------------------

# High-level orchestrator: reads FASTAs, validates lengths, builds position/state sets, and compares them.

class Alignments:

	# Initialize with verbosity flag and list of FASTA paths; internal containers are prepared lazily.

	def __init__(s, verbose: bool, fastas: List[str]):
		s.verbose = verbose
		s.fasta_paths = fastas

	# Read all FASTAs into dictionaries {seqID: SEQ}, handling .gz files and basic FASTA validation.

	def readFasta(s) -> None:
		s.fasta_dics: Dict[str, Dict[str, str]] = {}
		for file_name in s.fasta_paths:
			if s.verbose:
				pe(f"# Reading from {file_name}")
			seqs: Dict[str, str] = {}
			path = Path(file_name)
			try:
				if not path.exists():
					raise FileNotFoundError(f"File not found: {path}")
				if path.is_dir():
					raise IsADirectoryError(f"Is a directory, not a file: {path}")
				open_fn = gzip.open if path.suffix == ".gz" else open
				with open_fn(path, "r", encoding="utf-8", errors="strict") as handle:
					for rec in SeqIO.parse(handle, "fasta"):
						if rec.id in seqs:
							raise ValueError(f"Duplicate sequence ID '{rec.id}' in {path}")
						seqs[rec.id] = handleDashes(str(rec.seq)).upper()
				if not seqs:
					raise ValueError(f"No FASTA records found in {path}")
				s.fasta_dics[str(path)] = seqs
			except (FileNotFoundError, PermissionError, IsADirectoryError) as e:
				raise IOError(f"Cannot read '{path}': {e}") from e
			except UnicodeDecodeError as e:
				raise IOError(f"Encoding error while reading '{path}': {e}") from e
			except ValueError as e:
				raise IOError(f"Invalid FASTA content in '{path}': {e}") from e
			except Exception as e:
				raise IOError(f"Unexpected error while reading '{path}': {e.__class__.__name__}: {e}") from e
		return

	# Sanity-check that all alignments share identical sequence IDs and column counts.

	def checkAlignments(s) -> None:
		all_seqids: List[Set[str]] = []
		all_seqs: Dict[str, str] = {}
		for path in s.fasta_dics:
			alignment = s.fasta_dics[path]
			sizes: List[int] = []
			seqids: List[str] = []
			for seqid in alignment:
				seq = alignment[seqid]
				clean_seq = cleanSeq(seq)
				if not seqid in all_seqs:
					all_seqs[seqid] = clean_seq
				else:
					if all_seqs[seqid] != clean_seq:
						raise Exception(f"Sequences of terminal {seqid} are not the same in all alignments")
				sizes.append(len(seq))
				seqids.append(seqid)
			if len(set(sizes)) != 1:
				raise Exception(f"Sequences in {path} are of different lengths")
			all_seqids.append(set(seqids))
		if len({frozenset(s) for s in all_seqids}) != 1:
			raise Expection(f"The alignments do not share the same sequence identifiers")
		return
	
	# Build, for each alignment, the set-of-sets representing taxa sharing the same character at each site.

	def getPositions(s) -> None:
		fasta_positions_individuals: Dict[str, Dict[str, Set[str]]] = {}
		for path in s.fasta_dics:
			alignment = s.fasta_dics[path]
			fasta_positions_individuals[path]: Dict[str, Set[str]] = {}
			for seqid in alignment:
				seq = alignment[seqid]
				pos = seq2pos(seqid, seq)
				fasta_positions_individuals[path][seqid] = pos
		s.fasta_positions: Dict[str, Set[str]] = {}
		for path in fasta_positions_individuals:
			positions_basked: Dict[str, List[str]] = {}
			s.fasta_positions[path]: List[str] = []
			for seqid in fasta_positions_individuals[path]:
				positions = fasta_positions_individuals[path][seqid]
				count_position = 0
				for item in positions:
					if not count_position in positions_basked:
						positions_basked[count_position] = []
					positions_basked[count_position].append(item)
					count_position += 1
			for count_position in range(len(positions_basked)):
				s.fasta_positions[path].append(set(positions_basked[count_position]))
		return

	# Build, for each alignment, the set-of-sets representing taxa grouped by character state across sites.

	def getStates(s) -> None:
		fasta_states_individuals: Dict[str, Dict[str, Set[str]]] = {}
		for path in s.fasta_dics:
			alignment = s.fasta_dics[path]
			fasta_states_individuals[path]: Dict[str, Set[str]] = {}
			for seqid in alignment:
				seq = alignment[seqid]
				character_states = seq2char(seqid, seq)
				fasta_states_individuals[path][seqid] = character_states
		s.fasta_states: Dict[str, Set[str]] = {}
		for path in fasta_states_individuals:
			states_basked: Dict[str, List[str]] = {}
			s.fasta_states[path]: List[str] = []
			for seqid in fasta_states_individuals[path]:
				states = fasta_states_individuals[path][seqid]
				count_position = 0
				for item in states:
					if not count_position in states_basked:
						states_basked[count_position] = []
					states_basked[count_position].append(item)
					count_position += 1
			for count_position in range(len(states_basked)):
				s.fasta_states[path].append(set(states_basked[count_position]))
		return

	# Compare position-wise sets: report shared/unique counts and optimal matching similarity.

	def positionComparisons(s) -> None:
		po("File1,File2,Shared,Uniq1,Uniq2,PreliminaryScore,FinalScore")
		file_names = [fname for fname in s.fasta_positions]
		for fname1, fname2 in list(combinations(file_names, 2)):
			if s.verbose:
				pe(f"# Comparing alignment positions in {fname1} and {fname2}")
			positions1 = s.fasta_positions[fname1]
			positions2 = s.fasta_positions[fname2]
			set1 = {frozenset(s) for s in positions1}
			set2 = {frozenset(s) for s in positions2}
			shared = set1 & set2
			count_shared = len(shared)
			unique_to_positions1 = set1 - set2
			count_uniq1 = len(unique_to_positions1)
			unique_to_positions2 = set2 - set1
			count_uniq2 = len(unique_to_positions2)
			match_results = optimal_match_score([set(s) for s in unique_to_positions1], [set(s) for s in unique_to_positions2])
			score = match_results[0]
			final_score = final_similarity(count_shared, count_uniq1, count_uniq2, score)
			po(f"{fname1},{fname2},{count_shared},{count_uniq1},{count_uniq2},{score},{final_score}")
		return

	# Compare state-wise sets: report shared/unique counts and optimal matching similarity.

	def statesComparisons(s) -> None:
		po("File1,File2,Shared,Uniq1,Uniq2,PreliminaryScoreStates,FinalScoreStates")
		file_names = [fname for fname in s.fasta_states]
		for fname1, fname2 in list(combinations(file_names, 2)):
			if s.verbose:
				pe(f"# Comparing alignment states in {fname1} and {fname2}")
			states1 = s.fasta_states[fname1]
			states2 = s.fasta_states[fname2]
			set1 = {frozenset(s) for s in states1}
			set2 = {frozenset(s) for s in states2}
			shared = set1 & set2
			count_shared = len(shared)
			unique_to_states1 = set1 - set2
			count_uniq1 = len(unique_to_states1)
			unique_to_states2 = set2 - set1
			count_uniq2 = len(unique_to_states2)
			match_results = optimal_match_score([set(s) for s in unique_to_states1], [set(s) for s in unique_to_states2])
			score = match_results[0]
			final_score = final_similarity(count_shared, count_uniq1, count_uniq2, score)
			po(f"{fname1},{fname2},{count_shared},{count_uniq1},{count_uniq2},{score},{final_score}")
		return

# ---------------------------
# Main function
# ---------------------------

# Wire up the pipeline: read FASTAs → checks → build sets → run comparisons → print summaries.

def main(args: dict) -> None:
	a = Alignments(args.verbose, args.input)
	a.readFasta()
	a.checkAlignments()
	a.getPositions()
	a.positionComparisons()
	a.getStates()
	a.statesComparisons()
	return

# ---------------------------
# Run if executed as main
# ---------------------------

if __name__ == "__main__":
	args = get_args()
	main(args)