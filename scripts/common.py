"""Shared helpers for the FSJ recombination-mapping pipeline.

Ported from notebooks/scrubjay.ipynb so every pipeline stage (and any future
chromosome) uses the exact same logic the notebook validated on chromosome 1.
Nothing here is chromosome-specific: chromosome length is always read from
the tree sequence itself (`ts.sequence_length`), and chromosome identity is
passed in as a plain string ("1", "1A", "4A", "Z", ...) matching the FSJ
karyotype naming used in the .ped, GFF and linkage-map files.
"""

import csv
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import tskit
import yaml

MISSING_PARENT_CODE = "0"


# ---------------------------------------------------------------------------
# .ped parsing / pedigree construction (msprime PedigreeBuilder input)
# ---------------------------------------------------------------------------
def read_ped(path):
    """Read the .ped file using its header row to map column names to positions."""
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        header = [h.strip().strip("#").strip() for h in header]
        required = [
            "FamilyID", "IndividualID", "PaternalID", "MaternalID",
            "Sex", "CoverageClass", "Time", "IsFounder",
        ]
        missing = [c for c in required if c not in header]
        if missing:
            raise ValueError(f"Missing expected columns in .ped header: {missing}")
        idx = {name: header.index(name) for name in required}

        rows = []
        for line in reader:
            if not line or all(not v.strip() for v in line):
                continue
            rows.append({
                "fam": line[idx["FamilyID"]].strip(),
                "ind": line[idx["IndividualID"]].strip(),
                "pat": line[idx["PaternalID"]].strip(),
                "mat": line[idx["MaternalID"]].strip(),
                "sex": line[idx["Sex"]].strip(),
                "coverage_class": line[idx["CoverageClass"]].strip(),
                "time": float(line[idx["Time"]].strip()),
                "is_founder": line[idx["IsFounder"]].strip(),
            })
        return rows


def build_graph(ped_rows):
    """Build lookup dicts keyed by (fam, ind)."""

    def key(fam, ind):
        return (fam, ind)

    all_keys = set()
    row_of = {}
    parent_of = {}

    for row in ped_rows:
        k = key(row["fam"], row["ind"])
        all_keys.add(k)
        row_of[k] = row
        p0 = key(row["fam"], row["pat"]) if row["pat"] != MISSING_PARENT_CODE else None
        p1 = key(row["fam"], row["mat"]) if row["mat"] != MISSING_PARENT_CODE else None
        parent_of[k] = (p0, p1)

    return all_keys, row_of, parent_of


def sanity_check_founders(all_keys, row_of, parent_of):
    """Cross-check the IsFounder column against the parent columns."""
    problems = []
    for k in all_keys:
        p0, p1 = parent_of.get(k, (None, None))
        has_parents = p0 is not None or p1 is not None
        claimed_founder = row_of[k]["is_founder"].strip().lower() in ("1", "true", "yes", "t")
        if claimed_founder and has_parents:
            problems.append(f"{k}: marked IsFounder but has a parent listed")
        if not claimed_founder and not has_parents:
            problems.append(f"{k}: has no parents listed but IsFounder is false")
    if problems:
        print("Warning: IsFounder / parent-column mismatches found:", file=sys.stderr)
        for p in problems[:20]:
            print(f"  - {p}", file=sys.stderr)
        if len(problems) > 20:
            print(f"  ...and {len(problems) - 20} more", file=sys.stderr)


def validate_times(all_keys, row_of, parent_of):
    """msprime requires an individual's time to be strictly less than both parents'."""
    bad = []
    for k in all_keys:
        p0, p1 = parent_of.get(k, (None, None))
        t = row_of[k]["time"]
        for p in (p0, p1):
            if p is not None and row_of[p]["time"] <= t:
                bad.append(
                    f"{k} (time={t}) has parent {p} (time={row_of[p]['time']}), "
                    "parent time must be strictly greater"
                )
    if bad:
        raise ValueError(
            "Found individuals whose parent(s) are not strictly older:\n" + "\n".join(bad)
        )


def topological_order(all_keys, parent_of):
    """Sort individuals so every parent appears before its children (Kahn's algorithm)."""
    sorted_keys = sorted(all_keys)

    in_degree = {k: 0 for k in sorted_keys}
    for k in sorted_keys:
        p0, p1 = parent_of.get(k, (None, None))
        in_degree[k] = (1 if p0 else 0) + (1 if p1 else 0)

    children = {k: [] for k in sorted_keys}
    for k, (p0, p1) in parent_of.items():
        if p0:
            children[p0].append(k)
        if p1:
            children[p1].append(k)

    queue = sorted(k for k in sorted_keys if in_degree[k] == 0)
    order = []
    while queue:
        k = queue.pop(0)
        order.append(k)
        ready = []
        for c in children.get(k, []):
            in_degree[c] -= 1
            if in_degree[c] == 0:
                ready.append(c)
        queue.extend(sorted(ready))
        queue.sort()

    if len(order) != len(sorted_keys):
        raise ValueError("Cycle detected — cannot topologically sort pedigree")
    return order, children


def build_pedigree(ped_rows):
    import msprime

    all_keys, row_of, parent_of = build_graph(ped_rows)
    sanity_check_founders(all_keys, row_of, parent_of)
    validate_times(all_keys, row_of, parent_of)
    order, children = topological_order(all_keys, parent_of)

    pb = msprime.PedigreeBuilder(
        individuals_metadata_schema=tskit.MetadataSchema.permissive_json()
    )
    builder_id = {}

    for k in order:
        row = row_of[k]
        p0, p1 = parent_of.get(k, (None, None))
        parents = [builder_id[p0] if p0 else -1, builder_id[p1] if p1 else -1]
        is_sample = len(children.get(k, [])) == 0
        builder_id[k] = pb.add_individual(
            time=row["time"],
            parents=parents,
            is_sample=is_sample,
            metadata={"individual_id": row["ind"], "sex": row["sex"]},
        )

    return pb


def build_and_save_pedigree(ped_path, pedigree_out_path, sequence_length):
    """Parse the .ped file, build the pedigree, and save it so it never needs
    to be rebuilt from the raw .ped file again for this chromosome length."""
    ped_rows = read_ped(ped_path)
    pb = build_pedigree(ped_rows)
    pedigree_tables = pb.finalise(sequence_length=sequence_length)
    pedigree_tables.tree_sequence().dump(pedigree_out_path)
    print(f"Saved base pedigree ({pedigree_tables.individuals.num_rows} individuals) to {pedigree_out_path}")
    return pedigree_tables


# ---------------------------------------------------------------------------
# Individual ID <-> tree-sequence ID lookups
# ---------------------------------------------------------------------------
def get_id(ind):
    """FSJ pedigree ID for an individual, from either real (raw yaml-ish) or
    null/simulated (schema-decoded dict) metadata."""
    md = ind.metadata
    if isinstance(md, dict):
        return md.get("name") or md.get("individual_id") or md.get("id")
    s = md.decode("utf-8", "replace") if isinstance(md, (bytes, bytearray)) else str(md)
    for line in s.splitlines():
        if line.startswith("name:"):
            return line.split(":", 1)[1].strip()
        if line.startswith("individual_id:"):
            return line.split(":", 1)[1].strip()
    return None


def build_id_table(ts):
    """Lookup between .ped IndividualID and tree-sequence individual ID, plus
    per-individual info (sample status, founder status, time).

    Returns
    -------
    df : pandas.DataFrame   columns: ped_id, ts_id, is_sample, is_founder, time
    ped_to_ts : dict        {ped_individual_id: ts_individual_id}
    ts_to_ped : dict        {ts_individual_id: ped_individual_id}
    """
    rows = []
    ped_to_ts = {}
    ts_to_ped = {}

    for ind in ts.individuals():
        ped_id = get_id(ind)
        if ped_id is not None:
            ped_to_ts[ped_id] = ind.id
            ts_to_ped[ind.id] = ped_id

        is_sample = any(bool(ts.node(n).flags & tskit.NODE_IS_SAMPLE) for n in ind.nodes)
        is_founder = tuple(ind.parents) == (-1, -1)
        time = ts.node(ind.nodes[0]).time if len(ind.nodes) > 0 else None

        rows.append({
            "ped_id": ped_id, "ts_id": ind.id,
            "is_sample": is_sample, "is_founder": is_founder, "time": time,
        })

    df = pd.DataFrame(rows)
    return df, ped_to_ts, ts_to_ped


# ---------------------------------------------------------------------------
# Recombination-breakpoint counting
# ---------------------------------------------------------------------------
def per_node(ts):
    """Switches (S) and exposure (X) per node, used for CO/haplotype-ratio QC."""
    e = ts.tables.edges
    nt = ts.tables.nodes.time
    gap = nt[e.parent] - nt[e.child]
    expo = gap * (e.right - e.left) / 1e6
    order = np.lexsort((e.left, e.child))
    c, p = e.child[order], e.parent[order]
    sw = (c[1:] == c[:-1]) & (p[1:] != p[:-1])
    S = np.zeros(ts.num_nodes, np.int64)
    np.add.at(S, c[1:][sw], 1)
    X = np.zeros(ts.num_nodes)
    np.add.at(X, e.child, expo)
    return S, X


def co_per_gamete(ts):
    """Absolute crossovers per time-0 gamete — length-independent biological check."""
    nd, e = ts.tables.nodes, ts.tables.edges
    tips = np.array([n for n in ts.samples() if nd.time[n] == 0])
    co = np.array([(e.child == t).sum() - 1 for t in tips])
    return co, tips


def count_all_breakpoints(ts):
    """Breakpoint count per node, from a single pass over all trees."""
    n_nodes = ts.num_nodes
    breakpoint_counts = np.zeros(n_nodes, dtype=int)
    prev_parents = None

    for tree in ts.trees():
        parents = tree.parent_array[:n_nodes].copy()
        if prev_parents is not None:
            breakpoint_counts += parents != prev_parents
        prev_parents = parents

    return breakpoint_counts


def build_breakpoints_table_fast(ts, ped_to_ts):
    """Per-individual breakpoints table from a single precomputed pass over all trees."""
    node_breakpoints = count_all_breakpoints(ts)

    rows = []
    for ped_id, ts_id in ped_to_ts.items():
        individual = ts.individual(ts_id)
        node_ids = list(individual.nodes)
        row = {
            "ped_id": ped_id,
            "ts_id": ts_id,
            "total_breakpoints": sum(node_breakpoints[n] for n in node_ids),
        }
        for i, node_id in enumerate(node_ids):
            row[f"genome{i}_node_id"] = node_id
            row[f"genome{i}_breakpoints"] = int(node_breakpoints[node_id])
        rows.append(row)

    return pd.DataFrame(rows)


def informative_nodes(ts):
    """{node: (fsj_id, sorted breakpoint positions)} for every node with a parent edge."""
    e, nd = ts.tables.edges, ts.tables.nodes
    has_parent = np.zeros(ts.num_nodes, bool)
    has_parent[np.unique(e.child)] = True
    out = {}
    for ind in ts.individuals():
        fsj = get_id(ind)
        for node in ind.nodes:
            if has_parent[node]:
                b = np.sort(e.left[e.child == node])
                out[node] = (fsj, b[b > 0])
    return out


# ---------------------------------------------------------------------------
# Sibship / phasing-error filters
# ---------------------------------------------------------------------------
def read_sibships(ped_path, missing=MISSING_PARENT_CODE):
    """Return {individual: (father, mother)} and {(father,mother): [offspring]}
    for full-sib families (>=2 offspring)."""
    parents = {}
    with open(ped_path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        header = [h.strip().strip("#").strip() for h in next(r)]
        ix = {c: header.index(c) for c in ["IndividualID", "PaternalID", "MaternalID"]}
        for row in r:
            if not row or all(not v.strip() for v in row):
                continue
            parents[row[ix["IndividualID"]].strip()] = (
                row[ix["PaternalID"]].strip(), row[ix["MaternalID"]].strip())
    sibgroups = defaultdict(list)
    for ind, (pat, mat) in parents.items():
        if pat != missing and mat != missing:
            sibgroups[(pat, mat)].append(ind)
    sibgroups = {k: v for k, v in sibgroups.items() if len(v) >= 2}
    return parents, sibgroups


def _node_is_maternal(tseq, node_id):
    """True/False by parent-of-origin; None if the flag is unavailable for this node."""
    m = tseq.node(node_id).metadata
    if isinstance(m, dict):
        return m.get("is_maternal") != 0
    text = m.decode("utf-8", "replace") if isinstance(m, (bytes, bytearray)) else str(m)
    try:
        return yaml.safe_load(text).get("is_maternal", 0) != 0
    except Exception:
        return None


def is_maternal(tseq, node_id):
    meta = tseq.node(node_id).metadata
    if isinstance(meta, dict):
        return meta.get("is_maternal") != 0
    text = meta.decode("utf-8") if isinstance(meta, (bytes, bytearray)) else meta
    return yaml.safe_load(text)["is_maternal"] != 0


def sibling_recurrence_pairwise(nodes_by_id, sibgroups, tseq, window=1e4, role_by_node=None):
    """Flag breakpoints within `window` bp of a breakpoint on the SAME parental
    haplotype (paternal-vs-paternal, maternal-vs-maternal only) in a full sibling.
    Returns a long-format DataFrame: one row per (node, breakpoint), with
    n_sibs_matching, sibship, role."""
    rows = []
    for node, (fsj, bps) in nodes_by_id.items():
        if role_by_node is not None:
            role = role_by_node.get(node)
        else:
            is_mat = _node_is_maternal(tseq, node)
            role = "mat" if is_mat else ("pat" if is_mat is not None else None)
        for p in bps:
            rows.append({"fsj_id": fsj, "node": int(node), "position_bp": int(p), "role": role})
    d = pd.DataFrame(rows)
    d["position_Mb"] = (d.position_bp / 1e6).round(3)
    d["n_sibs_matching"] = 0
    d["sibship"] = None

    for (pat, mat), sibs in sibgroups.items():
        key = f"{pat}x{mat}"
        for role in ("pat", "mat"):
            role_mask = d.fsj_id.isin(sibs) & (d.role == role)
            pos_by_ind = {fsj: np.sort(grp["position_bp"].values)
                          for fsj, grp in d[role_mask].groupby("fsj_id")}
            present = [s for s in sibs if s in pos_by_ind]
            if len(present) < 2:
                continue
            d.loc[role_mask, "sibship"] = key
            for focal in present:
                focal_mask = role_mask & (d.fsj_id == focal)
                focal_pos = d.loc[focal_mask, "position_bp"].values
                match_counts = np.zeros(len(focal_pos), dtype=int)
                for other in present:
                    if other == focal:
                        continue
                    op = pos_by_ind[other]
                    idx = np.searchsorted(op, focal_pos)
                    left = np.clip(idx - 1, 0, len(op) - 1)
                    right = np.clip(idx, 0, len(op) - 1)
                    nearest = np.minimum(np.abs(op[left] - focal_pos), np.abs(op[right] - focal_pos))
                    match_counts += (nearest <= window).astype(int)
                d.loc[focal_mask, "n_sibs_matching"] = match_counts
    return d


def all_parental_switches(tseq, sibgroups, parents, id_to_ind):
    """{offspring: {'pat': {pos:dir}, 'mat': {pos:dir}}} — haplotype-switch
    direction per offspring, used to detect opposite-direction phasing errors."""
    e = tseq.tables.edges
    nd = tseq.tables.nodes
    parent_nodes_of = defaultdict(set)
    for c, p in zip(e.child, e.parent):
        parent_nodes_of[c].add(p)

    switches = defaultdict(lambda: {"pat": {}, "mat": {}})
    off_to_sib = {}
    for (pat, mat), sibs in sibgroups.items():
        for s in sibs:
            if s not in id_to_ind:
                continue
            off_to_sib[s] = f"{pat}x{mat}"
            for node in id_to_ind[s].nodes:
                pnodes = parent_nodes_of.get(node, set())
                par_names = {get_id(tseq.individual(nd.individual[pn]))
                             for pn in pnodes if nd.individual[pn] != -1}
                which = "pat" if pat in par_names else ("mat" if mat in par_names else None)
                if which is None:
                    continue
                haps = sorted(pnodes)
                if len(haps) < 2:
                    continue
                h0, h1 = haps[0], haps[1]
                mask = e.child == node
                lefts = e.left[mask]
                pars = e.parent[mask]
                order = np.argsort(lefts)
                lefts, pars = lefts[order], pars[order]
                sw = {}
                prev = None
                for L, p in zip(lefts, pars):
                    cur = 0 if p == h0 else (1 if p == h1 else None)
                    if prev is not None and cur is not None and cur != prev:
                        sw[int(L)] = cur - prev
                    if cur is not None:
                        prev = cur
                switches[s][which] = sw
    return switches, off_to_sib


def opposite_direction_flags(switches, off_to_sib, sibgroups, id_to_ind, window=1e4):
    """Flag breakpoints where one sibling switches P/M0 -> P/M1 and another
    switches P/M1 -> P/M0 at the same spot (a hallmark of a phasing error)."""
    results = []
    for (pat, mat), sibs in sibgroups.items():
        present = [s for s in sibs if s in switches]
        for i in range(len(present)):
            for j in range(i + 1, len(present)):
                a, b = present[i], present[j]
                for role in ("pat", "mat"):
                    sw_a = switches[a].get(role, {})
                    sw_b = switches[b].get(role, {})
                    if not sw_a or not sw_b:
                        continue
                    bpos = np.array(sorted(sw_b.keys()))
                    for pa, da in sw_a.items():
                        idx = np.searchsorted(bpos, pa)
                        for k in (idx - 1, idx):
                            if 0 <= k < len(bpos) and abs(bpos[k] - pa) <= window:
                                db = sw_b[int(bpos[k])]
                                if da == -db:
                                    results.append({
                                        "sibship": f"{pat}x{mat}", "parent": role,
                                        "sib_a": a, "pos_a": pa, "dir_a": da,
                                        "sib_b": b, "pos_b": int(bpos[k]), "dir_b": db})
    return pd.DataFrame(results)


def per_parent_recurrence(bps_by_node_with_id, parents, window=1e3, min_share=2):
    """Flag breakpoints recurring across offspring of the same FATHER or same
    MOTHER (catches half-sib recurrence that full-sib filtering misses)."""
    by_father, by_mother = defaultdict(list), defaultdict(list)
    for node, (fsj, bps) in bps_by_node_with_id.items():
        if fsj not in parents:
            continue
        pat, mat = parents[fsj]
        if pat != MISSING_PARENT_CODE:
            by_father[pat].append((fsj, bps))
        if mat != MISSING_PARENT_CODE:
            by_mother[mat].append((fsj, bps))

    flagged = defaultdict(set)
    for parent_group in (by_father, by_mother):
        for parent, offspring in parent_group.items():
            if len(offspring) < min_share:
                continue
            all_pos = []
            for fsj, bps in offspring:
                all_pos.extend((fsj, p) for p in bps)
            pos_arr = sorted(set(p for _, p in all_pos))
            for center in pos_arr:
                sharers = set(fsj for fsj, p in all_pos if abs(p - center) <= window)
                if len(sharers) >= min_share:
                    for fsj, p in all_pos:
                        if abs(p - center) <= window:
                            flagged[fsj].add(p)
    return flagged


# ---------------------------------------------------------------------------
# Distance-based filter (nearest-neighbour antimode)
# ---------------------------------------------------------------------------
def breakpoint_distance_table(bps_by_node, nodes_meta, chrom_length):
    """One row per breakpoint: individual, node, position, distance to nearest
    neighbour on the same haplotype (and to each side)."""
    rows = []
    for node, bps in bps_by_node.items():
        fsj = nodes_meta[node][0]
        bps = np.sort(np.unique(bps))
        if len(bps) == 0:
            continue
        left = np.empty(len(bps))
        right = np.empty(len(bps))
        left[0] = bps[0]
        left[1:] = np.diff(bps)
        right[-1] = chrom_length - bps[-1]
        right[:-1] = np.diff(bps)
        nearest = np.minimum(left, right)

        for p, l, r, n in zip(bps, left, right, nearest):
            rows.append({
                "fsj_id": fsj, "node": int(node), "position_Mb": round(p / 1e6, 2),
                "dist_left_bp": round(l), "dist_right_bp": r, "nearest_bp": n,
            })
    df = pd.DataFrame(rows).sort_values(["fsj_id", "node"])
    return df.reset_index(drop=True)


def find_antimode_cutoff(nearest_bp, bins=80, smooth_window=7, peak_floor_frac=0.02):
    """Find the antimode (valley) of the bimodal nearest-neighbour-distance
    distribution, separating "true" spacing from clustered/artifactual
    breakpoints. Returns the cutoff in bp, or None if no clear bimodal
    structure is found (falls back to no distance filtering upstream)."""
    d = nearest_bp[nearest_bp > 0]
    if len(d) == 0:
        return None
    logd = np.log10(d)
    counts, edges = np.histogram(logd, bins=bins, density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    smooth = np.convolve(counts, np.ones(smooth_window) / smooth_window, mode="same")
    peaks = np.where((smooth[1:-1] > smooth[:-2]) & (smooth[1:-1] > smooth[2:]))[0] + 1
    prom = [p for p in peaks if smooth[p] > peak_floor_frac * smooth.max()]
    if len(prom) < 2:
        return None
    lo, hi = prom[0], prom[-1]
    valley = lo + np.argmin(smooth[lo:hi])
    return float(10 ** centers[valley])


def nearest_neighbor_filter(bps_by_node, min_dist, chrom_length):
    kept = {}
    for node, bps in bps_by_node.items():
        bps = np.sort(np.unique(bps))
        if len(bps) < 2:
            kept[node] = bps
            continue
        left = np.empty(len(bps))
        right = np.empty(len(bps))
        left[0] = bps[0]
        left[1:] = np.diff(bps)
        right[-1] = chrom_length - bps[-1]
        right[:-1] = np.diff(bps)
        nearest = np.minimum(left, right)
        kept[node] = bps[nearest >= min_dist]
    return kept


# ---------------------------------------------------------------------------
# Marey maps
# ---------------------------------------------------------------------------
def build_mareymap(bps_by_node_sexfilter, chrom_len, n_bins=500):
    """Cumulative cM vs bp from per-gamete crossover positions. Genetic distance
    accrues as crossover density integrated along the chromosome."""
    edges = np.linspace(0, chrom_len, n_bins + 1)
    n_gametes = len(bps_by_node_sexfilter)
    all_bp = (np.concatenate([np.asarray(b) for b in bps_by_node_sexfilter.values()])
              if n_gametes else np.array([]))
    counts, _ = np.histogram(all_bp, bins=edges)
    rec_freq = counts / n_gametes if n_gametes else counts.astype(float)
    cum_cM = np.concatenate([[0], np.cumsum(rec_freq) * 100])
    return edges, cum_cM


# ---------------------------------------------------------------------------
# GFF annotation
# ---------------------------------------------------------------------------
def load_genes(gff_path, chrom):
    genes = []
    with open(gff_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            if c[0] != chrom:
                continue
            genes.append({"start": int(c[3]), "end": int(c[4]), "strand": c[6], "attr": c[8].strip()})
    return pd.DataFrame(genes).sort_values("start").reset_index(drop=True)


def load_features(gff_path, chrom, feature_types):
    """Load intervals for the given feature types (e.g. exon, CDS, gene, intron)."""
    feats = {ft: [] for ft in feature_types}
    with open(gff_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or "\t" not in line:
                continue
            c = line.split("\t")
            if len(c) < 9 or c[0] != chrom:
                continue
            if c[2] in feats:
                feats[c[2]].append((int(c[3]), int(c[4])))
    return feats


def merge_intervals(intervals):
    """Merge overlapping intervals, return sorted list and total bp covered."""
    if not intervals:
        return [], 0
    iv = sorted(intervals)
    merged = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    total = sum(e - s for s, e in merged)
    return merged, total


def frac_in(positions, merged):
    """Fraction of positions falling in a sorted merged-interval list (binary search)."""
    if not merged or len(positions) == 0:
        return 0, 0
    starts = np.array([s for s, e in merged])
    ends = np.array([e for s, e in merged])
    hit = 0
    for p in positions:
        i = np.searchsorted(starts, p) - 1
        if 0 <= i < len(ends) and p <= ends[i]:
            hit += 1
    return hit, hit / len(positions)


def dist_to(pos, sites):
    """Distance from each position to the nearest site in a sorted array."""
    assert np.all(np.diff(sites) >= 0), "sites must be sorted"
    i = np.searchsorted(sites, pos)
    lo = sites[np.clip(i - 1, 0, len(sites) - 1)]
    hi = sites[np.clip(i, 0, len(sites) - 1)]
    return np.minimum(np.abs(pos - lo), np.abs(pos - hi))


def pooled_breakpoints(bps_by_node):
    """Flatten a {node: positions} or {node: (fsj_id, positions)} dict (as
    returned by informative_nodes()) into one sorted array."""
    arrs = []
    for v in bps_by_node.values():
        b = v[1] if isinstance(v, tuple) else v
        b = np.asarray(b)
        if b.size:
            arrs.append(b)
    return np.sort(np.concatenate(arrs)).astype(np.int64) if arrs else np.array([], np.int64)
