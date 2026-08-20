"""Standalone IBD-tract plotter for a PedigreeTS (Widget 2 of the pedigree viewer).

Self-contained: needs only `tskit` + `matplotlib`. Given a PedigreeTS and a focal
individual ID, it simplifies down to the six trio nodes (focal chrom1/2, father 0/1,
mother 0/1) and paints, per local tree, which parental haplotype each of the focal's
two chromosomes matches — dark/light blue = father's two haplotypes (P0/P1),
dark/light orange = mother's two haplotypes (M0/M1).

The trio is derived straight from the tree sequence:
  * individuals are keyed by their metadata `name` / `file_id`;
  * the focal's paternal vs maternal haplotype nodes are read from the node
    `is_maternal` flag;
  * father/mother are identified as the individuals owning the tree-parent of the
    focal's paternal / maternal node (so no external pedigree file is required).

Usage
-----
    import tskit
    from plot_ibd_tracts import plot_ibd

    ts = tskit.load("1.shapeit5.hmm.v3.pedigree.trees")
    plot_ibd(ts, "FSJ2252")                       # show interactively
    plot_ibd(ts, "FSJ2252", out_path="ibd.png")   # save to file
    plot_ibd(ts, "FSJ2252", region=(2_000_000, 8_000_000))  # zoom (bp)

Or from the command line:
    python plot_ibd_tracts.py <trees_file> <focal_id> [out.png]
"""

import tskit
import yaml
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# ---------------------------------------------------------------------------
# Metadata helpers (ported from projectutils / app_core)
# ---------------------------------------------------------------------------
def get_name_from_metadata(metadata):
    """Extract the individual name from individual/node metadata.

    Real (phased) tree sequences carry `file_id`; simulated/null tree
    sequences built from the .ped pedigree carry `individual_id` instead.
    """
    if isinstance(metadata, dict):
        return metadata.get("file_id") or metadata.get("individual_id") or metadata.get("name")
    elif isinstance(metadata, (bytes, bytearray)):
        return yaml.safe_load(metadata.decode("utf-8"))["name"]
    elif isinstance(metadata, str):
        return yaml.safe_load(metadata)["name"]
    else:
        raise ValueError("Unknown metadata type: {}".format(type(metadata)))


def _is_maternal(ts, node_id):
    """True if `node_id` is its individual's maternal haplotype.

    Real tree sequences flag this directly in node metadata (`is_maternal`).
    Simulated/null tree sequences carry no node metadata at all, so fall back
    to the individual's pedigree `parents` ([father, mother]) plus tree
    topology: whichever pedigree parent owns this node's tree-parent tells us
    which haplotype it is.
    """
    meta = ts.node(node_id).metadata
    if isinstance(meta, dict) and "is_maternal" in meta:
        return meta["is_maternal"] != 0
    if isinstance(meta, (bytes, bytearray, str)) and meta:
        text = meta.decode("utf-8") if isinstance(meta, (bytes, bytearray)) else meta
        parsed = yaml.safe_load(text)
        if parsed and "is_maternal" in parsed:
            return parsed["is_maternal"] != 0

    ind = ts.individual(ts.node(node_id).individual)
    if len(ind.parents) != 2:
        raise ValueError(
            f"node {node_id}: no is_maternal metadata, and its individual has "
            f"no pedigree parents to fall back on"
        )
    father, mother = ind.parents
    parent_ind = _individual_of(ts, node_id)
    if parent_ind == mother:
        return True
    if parent_ind == father:
        return False
    raise ValueError(
        f"node {node_id}: tree-parent individual {parent_ind} matches neither "
        f"pedigree parent (father={father}, mother={mother})"
    )


# ---------------------------------------------------------------------------
# Derive the six trio nodes from the tree sequence
# ---------------------------------------------------------------------------
def _individual_of(ts, node_id):
    """The individual owning the tree-parent of `node_id` (scans trees until found)."""
    for tree in ts.trees():
        p = tree.parent(node_id)
        if p != tskit.NULL and ts.node(p).individual != tskit.NULL:
            return ts.node(p).individual
    return tskit.NULL


def derive_trio_nodes(ts, focal_id):
    """Return the six trio node IDs [focal_pat, focal_mat, dad0, dad1, mum0, mum1].

    Father / mother are inferred from the tree topology: the individual owning the
    tree-parent of the focal's paternal (resp. maternal) haplotype node.
    """
    name_to_ind = {get_name_from_metadata(i.metadata): i for i in ts.individuals()}
    if focal_id not in name_to_ind:
        raise KeyError(f"{focal_id!r} not found in the tree sequence")
    focal = name_to_ind[focal_id]

    pat_node = mat_node = None
    for n in focal.nodes:
        if _is_maternal(ts, n):
            mat_node = n
        else:
            pat_node = n
    if pat_node is None or mat_node is None:
        raise KeyError(f"{focal_id} is missing a paternal/maternal node in the TS")

    dad = _individual_of(ts, pat_node)
    mum = _individual_of(ts, mat_node)
    if dad == tskit.NULL or mum == tskit.NULL:
        raise ValueError(
            f"{focal_id} has no parental haplotype nodes in the TS "
            "(is it a founder?) — cannot build a trio IBD plot"
        )

    six = (
        list(focal.nodes)
        + list(ts.individual(dad).nodes)
        + list(ts.individual(mum).nodes)
    )
    return six


# ---------------------------------------------------------------------------
# IBD rendering (ported from tsinferutils.SampleIBDPlot, PedigreeTS variant)
# ---------------------------------------------------------------------------
# node 0/1 = focal chrom1/chrom2; 2/3 = father hap0/1; 4/5 = mother hap0/1
_NODE_COLOURS = {2: "C0", 3: "#66b3ff", 4: "C1", 5: "#ffcc99"}


def _nearest_neighbours(tree):
    """For each focal chromosome (nodes 0, 1), the parental node(s) it matches."""
    result = []
    for u in (0, 1):
        a = tree.parent(u)
        if a in {2, 3, 4, 5}:
            result.append([a])
        elif a != tskit.NULL:
            result.append([c for c in tree.leaves(a) if c > 1])
        else:
            result.append([])
    return result

def _tracts(trio):
    """Per focal haplotype: [[left, right, neighbours_tuple], ...], runs coalesced."""
    out = [[], []]
    for t in trio.trees():
        nn = _nearest_neighbours(t)
        l, r = t.interval
        for h in (0, 1):
            key = tuple(sorted(n for n in nn[h] if n in _NODE_COLOURS))
            if out[h] and out[h][-1][2] == key:
                out[h][-1][1] = r
            else:
                out[h].append([l, r, key])
    return out


def _merge_excluded(tracts, positions, tol=1.0):
    """Drop tract boundaries within tol bp of a flagged position; longer flank wins."""
    if not positions or not tracts:
        return tracts, 0
    merged, removed = [list(tracts[0])], 0
    for l, r, key in tracts[1:]:
        b = merged[-1][1]
        if any(abs(b - p) <= tol for p in positions):
            if (r - l) > (merged[-1][1] - merged[-1][0]):
                merged[-1][2] = key          # longer side keeps its label
            merged[-1][1] = r
            removed += 1
        else:
            merged.append([l, r, key])
    return merged, removed

def _paint(ax, loc_start, loc_end, seq_len, neighbours):
    neighbours = sorted(n for n in neighbours if n in _NODE_COLOURS)
    if len(neighbours) == 1:
        ax.add_patch(
            Rectangle(
                (loc_start / seq_len, 0),
                (loc_end - loc_start) / seq_len,
                1,
                color=_NODE_COLOURS[neighbours[0]],
            )
        )
    elif len(neighbours) == 2:
        h = 0.5
        for i, n in enumerate(neighbours):
            ax.add_patch(
                Rectangle(
                    (loc_start / seq_len, i * h),
                    (loc_end - loc_start) / seq_len,
                    h,
                    color=_NODE_COLOURS[n],
                )
            )


def plot_ibd(
    ts, focal_id, out_path=None, region=None, title=None,
    chrom=1, centromere=None,
    exclude_positions=None,      # iterable of bp to merge away
    exclude_tol=1.0,             # match tolerance in bp
):
    """Plot PedigreeTS IBD tracts for `focal_id`.

    Parameters
    ----------
    ts : tskit.TreeSequence | str
        A loaded PedigreeTS or a path to a `.trees` file.
    focal_id : str
        Individual name/token (e.g. "FSJ2252").
    out_path : str | None
        If given, save a PNG here; otherwise show interactively.
    region : (int, int) | None
        Zoom the x-axis to [start, end] in bp.
    title : str | None
        Override the plot title.
    chrom : int | str
        Chromosome label used in the default title.
    centromere : (int, int) | None
        If given, shade this bp interval grey on both chromosomes.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(ts, str):
        ts = tskit.load(ts)

    trio = ts.simplify(derive_trio_nodes(ts, focal_id))
    sl = trio.sequence_length
    r_start, r_end = region if region else (None, None)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 4), sharex=True)
    plt.subplots_adjust(hspace=0.3)
    for ax, label in ((ax1, "Chrom1"), (ax2, "Chrom2")):
        ax.get_yaxis().set_visible(False)
        ax.text(
            -0.02, 0.5, label, fontsize=12, ha="center", va="center",
            rotation="vertical", transform=ax.transAxes,
        )

    tr1, tr2 = _tracts(trio)
    excl = set(exclude_positions or ())
    tr1, rm1 = _merge_excluded(tr1, excl, exclude_tol)
    tr2, rm2 = _merge_excluded(tr2, excl, exclude_tol)
    if excl:
        print(f"{focal_id}: {len(excl)} flagged positions -> "
              f"removed {rm1} paternal + {rm2} maternal crossovers")
    for l, r, key in tr1:
        _paint(ax1, l, r, sl, key)
    for l, r, key in tr2:
        _paint(ax2, l, r, sl, key)

    # Legend to the right of the plot.
    lax = fig.add_axes([0.95, 0.1, 0.05, 0.8])
    lax.axis("off")
    lax.legend(
        handles=[Rectangle((0, 0), 1, 1, color=_NODE_COLOURS[n]) for n in (2, 3, 4, 5)],
        labels=["P0", "P1", "M0", "M1"],
        fontsize=10,
    )

    if centromere is not None:
        c0, c1 = centromere
        for ax in (ax1, ax2):
            ax.add_patch(
                Rectangle((c0 / sl, 0), (c1 - c0) / sl, 1, color="grey", alpha=0.8)
            )

    zoomed = r_start is not None and r_end is not None
    if zoomed:
        ax1.set_xlim(r_start / sl, r_end / sl)
        span = r_end - r_start
        for step in (5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7):
            if span / step <= 10:
                break
        step = int(step)
        ticks = range((int(r_start) // step + 1) * step, int(r_end), step)
        ax2.set_xticks([bp / sl for bp in ticks])
        ax2.set_xticklabels([f"{bp / 1e6:.1f} Mb" for bp in ticks], fontsize=10)
    else:
        step = 10_000_000
        ticks = range(step, int(sl), step)
        ax2.set_xticks([bp / sl for bp in ticks])
        ax2.set_xticklabels([f"{bp / 1e6:.0f}" for bp in ticks], fontsize=10)
        ax2.set_xlabel("Position (Mb)", fontsize=12)

    if title is None:
        zoom_str = f"  [{r_start:,}-{r_end:,} bp]" if zoomed else ""
        title = f"PedigreeTS IBD - {focal_id}  chr{chrom}{zoom_str}"
    plt.suptitle(title, fontsize=16, fontweight="bold")

    if out_path is not None:
        plt.savefig(out_path)
    else:
        plt.show()
    return fig


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        sys.exit("usage: python plot_ibd_tracts.py <trees_file> <focal_id> [out.png]")
    trees_file, focal = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else None
    plot_ibd(trees_file, focal, out_path=out)
    if out:
        print(f"wrote {out}")