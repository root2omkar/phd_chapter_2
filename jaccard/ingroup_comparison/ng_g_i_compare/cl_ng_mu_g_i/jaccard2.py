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

def get_args():
    """
    Parse command-line arguments:
    -i / --input: list of FASTA files
    -v / --verbose: enable verbose output
    """
    parser = argparse.ArgumentParser(
        description="Reads different alignments in FASTA format and compares them."
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        nargs='+',
        help="List FASTA files separated by space",
        required=True
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Increase output verbosity"
    )
    return parser.parse_args()

def po(message: str) -> None:
    """Print a message to stdout (normal output)."""
    sys.stdout.write(f"{message}\n")

def pe(message: str) -> None:
    """Print a message to stderr (warnings/errors/verbose)."""
    sys.stderr.write(f"{message}\n")

def handleDashes(seq: str) -> str:
    """
    Normalize gaps in sequences:
    - Replace all dash variants with '-'
    - Replace terminal Ns with '-' to treat them as gaps
    """
    new_seq = re.sub(r"[-–—]", "-", seq)
    new_seq = re.sub(r'^N+', lambda m: '-' * len(m.group(0)), new_seq)
    new_seq = re.sub(r'N+$', lambda m: '-' * len(m.group(0)), new_seq)
    return new_seq

def seq2pos(seqid: str, seq: str) -> List[str]:
    """
    Encode each sequence position as a unique key:
    - Gaps are labeled as 'seqID:-'
    - Non-gaps labeled as 'seqID:position_index'
    Used to build position-wise sets for Jaccard comparisons.
    """
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

def seq2char(seqid: str, seq: str) -> List[str]:
    """
    Encode each character state with sequence ID:
    'seqID:char'
    Used to build state-wise sets for Jaccard comparisons.
    """
    list_states: List[str] = []
    for char in seq:
        list_states.append(f"{seqid}:{char}")
    return list_states

def jaccard(a: Set[str], b: Set[str]) -> float:
    """
    Compute Jaccard similarity between two sets:
    |A ∩ B| / |A ∪ B|
    Returns 1.0 if both sets are empty.
    """
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0

def hungarian(cost):
    """
    Hungarian algorithm (pure Python implementation)
    Solves minimal cost assignment problem (square cost matrix required).
    Returns list of assigned columns for each row.
    """
    n = len(cost)
    u = [0.0] * (n + 1)
    v = [0.0] * (n + 1)
    p = [0] * (n + 1)
    way = [0] * (n + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [math.inf] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = math.inf
            j1 = 0
            for j in range(1, n + 1):
                if not used[j]:
                    cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
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
    match = [0] * (n + 1)
    for j in range(1, n + 1):
        match[p[j]] = j
    return match[1:]

def optimal_match_score(U: List[Set[str]], V: List[Set[str]]) -> Tuple[float, list]:
    """
    Compute optimal one-to-one matching between two lists of sets,
    maximizing total Jaccard similarity. Returns (score, matching).
    """
    m, n = len(U), len(V)
    k = max(m, n)
    S = [[0.0] * k for _ in range(k)]
    for i in range(m):
        for j in range(n):
            S[i][j] = jaccard(U[i], V[j])
    C = [[1.0 - S[i][j] for j in range(k)] for i in range(k)]
    assignment = hungarian(C)
    pairs = []
    sim_sum = 0.0
    for i in range(k):
        j = assignment[i] - 1
        if i < m and j < n:
            sim = S[i][j]
            pairs.append((i, j, sim))
            sim_sum += sim
        else:
            pairs.append((i, None if j >= n else j, 0.0))
    score = sim_sum / k if k else 1.0
    return score, pairs

def final_similarity(shared_count: int, m: int, n: int, match_score: float) -> float:
    """
    Combine exact matches and optimal matching score into a single similarity [0,1].
    """
    k = max(m, n)
    total = shared_count + m + n
    if total == 0:
        return 1.0
    return (float(shared_count) + match_score * float(k)) / float(total)

# ---------------------------
# Alignments Class
# ---------------------------

class Alignments:

    def __init__(s, verbose: bool, fastas: List[str]):
        s.verbose = verbose
        s.fasta_paths = fastas

    def readFasta(s) -> None:
        """
        Read all FASTA files into dictionaries {seqID: sequence}.
        Handles gzip files and validates basic FASTA structure.
        """
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
            except Exception as e:
                raise IOError(f"Error reading '{path}': {e}") from e

    def checkAlignments(s) -> None:
        """
        Ensure all alignments:
        1. Have the same sequence IDs.
        2. Have sequences of the same length in each alignment.
        """
        all_seqids: List[Set[str]] = []

        for path, alignment in s.fasta_dics.items():
            sizes = [len(seq) for seq in alignment.values()]
            seqids = set(alignment.keys())

            # Check that all sequences in this alignment have the same length
            if len(set(sizes)) != 1:
                raise Exception(f"Sequences in {path} have different lengths")

            all_seqids.append(seqids)

        # Check that all alignments have the same set of sequence IDs
        if len({frozenset(s) for s in all_seqids}) != 1:
            raise Exception("The alignments do not share the same sequence identifiers")

    def getPositions(s) -> None:
        """
        Build position-wise sets for each alignment.
        Each position is a set of sequence identifiers sharing that position.
        """
        s.fasta_positions: Dict[str, List[Set[str]]] = {}
        for path, alignment in s.fasta_dics.items():
            positions_basked: Dict[int, List[str]] = {}
            for seqid, seq in alignment.items():
                pos_list = seq2pos(seqid, seq)
                for idx, code in enumerate(pos_list):
                    positions_basked.setdefault(idx, []).append(code)
            s.fasta_positions[path] = [set(positions_basked[i]) for i in range(len(positions_basked))]

    def getStates(s) -> None:
        """
        Build state-wise sets for each alignment.
        Each set groups taxa by character state at a site.
        """
        s.fasta_states: Dict[str, List[Set[str]]] = {}
        for path, alignment in s.fasta_dics.items():
            states_basked: Dict[int, List[str]] = {}
            for seqid, seq in alignment.items():
                state_list = seq2char(seqid, seq)
                for idx, code in enumerate(state_list):
                    states_basked.setdefault(idx, []).append(code)
            s.fasta_states[path] = [set(states_basked[i]) for i in range(len(states_basked))]

    def positionComparisons(s) -> None:
        """
        Compare position-wise sets between alignments using Jaccard similarity.
        Reports shared/unique counts and combined score.
        """
        po("File1,File2,Shared,Uniq1,Uniq2,PreliminaryScore,FinalScore")
        files = list(s.fasta_positions.keys())
        for f1, f2 in combinations(files, 2):
            positions1 = s.fasta_positions[f1]
            positions2 = s.fasta_positions[f2]
            set1 = {frozenset(p) for p in positions1}
            set2 = {frozenset(p) for p in positions2}
            shared = set1 & set2
            uniq1 = set1 - set2
            uniq2 = set2 - set1
            match_score, _ = optimal_match_score([set(s) for s in uniq1], [set(s) for s in uniq2])
            final_score = final_similarity(len(shared), len(uniq1), len(uniq2), match_score)
            po(f"{f1},{f2},{len(shared)},{len(uniq1)},{len(uniq2)},{match_score},{final_score}")

    def statesComparisons(s) -> None:
        """
        Compare state-wise sets between alignments using Jaccard similarity.
        Reports shared/unique counts and combined score.
        """
        po("File1,File2,Shared,Uniq1,Uniq2,PreliminaryScoreStates,FinalScoreStates")
        files = list(s.fasta_states.keys())
        for f1, f2 in combinations(files, 2):
            states1 = s.fasta_states[f1]
            states2 = s.fasta_states[f2]
            set1 = {frozenset(p) for p in states1}
            set2 = {frozenset(p) for p in states2}
            shared = set1 & set2
            uniq1 = set1 - set2
            uniq2 = set2 - set1
            match_score, _ = optimal_match_score([set(s) for s in uniq1], [set(s) for s in uniq2])
            final_score = final_similarity(len(shared), len(uniq1), len(uniq2), match_score)
            po(f"{f1},{f2},{len(shared)},{len(uniq1)},{len(uniq2)},{match_score},{final_score}")

# ---------------------------
# Main
# ---------------------------

def main(args) -> None:
    a = Alignments(args.verbose, args.input)
    a.readFasta()
    a.checkAlignments()
    a.getPositions()
    a.positionComparisons()
    a.getStates()
    a.statesComparisons()

if __name__ == "__main__":
    args = get_args()
    main(args)
