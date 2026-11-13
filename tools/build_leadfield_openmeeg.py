"""EEG leadfield of any electrode set on the ICBM152 template head, with OpenMEEG.

Same head model as the pipeline's FieldTrip route (get_ftHeadmodel_and_ftLeadfield):
three nested BEM surfaces from Brainstorm's ICBM152 anatomy (scalp, outer
skull, inner skull), conductivities brain/skull/scalp = 1 / 0.0125 / 1,
symmetric BEM; sources = the 15,002 vertices of tess_cortex_pial_low.mat
(dipoles along the cortical normal). Used when FieldTrip cannot run its
'openmeeg' method because the OpenMEEG command-line programs are missing:
the Python package of OpenMEEG gives the same gain (checked against
FieldTrip's leadfield: median r = 1.000000 column by column).

Input: a Brainstorm channel_*.mat or a CSV with columns label,x,y,z (metres,
Brainstorm ICBM152 coordinates). Electrodes must lie on the template scalp:
the script stops (check_positions) when the coordinates look like mm or cm,
when the electrodes are more than 1 cm from the scalp surface (median), or
when Fpz/T7/Oz/T8 are more than 30 degrees from their expected azimuth (a
rotated cap).
--force turns these errors into warnings. A cap without these labels (e.g.
ANT's equidistant Z1, LL1, ...) cannot be checked for rotation (a message
says so). After the BEM, a physics check (sign and localisation of the gain)
at Cz, Fpz, Oz, T7, T8, or, without these labels, at the electrodes nearest
to their template positions (within 2 cm).

Another head (S.anatomy, folders made by tools/make_anatomy.m): --anatomy
FOLDER takes tess_head.mat, tess_outerskull.mat and tess_innerskull.mat from
that folder and --cortex FILE the sources (default: FOLDER/tess_cortex.mat);
the electrodes are then those of FOLDER's channel file, on that head's scalp
(same checks). Without --anatomy: the pipeline's ICBM152 of inputs/.
OpenMEEG gives non-finite values for some exactly degenerate meshes (e.g.
Brainstorm's 1922-vertex ICBM152 skulls, a face with an exactly coplanar
neighbour): the gain is then computed again with the vertices as Brainstorm
writes them for OpenMEEG (out_tess_tri: millimetres, %g), and the output's
bem_vertices says so.

Output (--out): G_normal (nChannels x nSources, single, no reference), labels,
elec_pos, src_pos, src_nrm, provenance (anatomy: the surfaces used,
bem_vertices: their precision); with
--free also G_free (nChannels x 3*nSources, x/y/z interleaved, double), the
free-orientation gain.

Usage (from the pipeline folder; Python 3.10+, pip install openmeeg scipy numpy):
  python tools/build_leadfield_openmeeg.py --channel-file my_cap.csv --out "1- EEG source connectivity/inputs/ftLeadfield_my_cap_normal.mat"
"""
import argparse
import csv
import datetime
import os
import sys

import numpy as np
import openmeeg as om
import scipy.io as sio
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANAT = os.path.join(ROOT, "1- EEG source connectivity", "inputs")
CONDUCTIVITY = (1.0, 0.0125, 1.0)  # brain, skull, scalp (inner to outer)


def load_surface(name, folder=ANAT):
    m = sio.loadmat(os.path.join(folder, name), squeeze_me=True, struct_as_record=False)
    return np.asarray(m["Vertices"], dtype=float), np.asarray(m["Faces"], dtype=np.int64) - 1, m


def signed_volume(V, F):
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0


def orient_for_openmeeg(V, F, name):
    # With OpenMEEG 2.6's make_nested_geometry the triangles must have a negative
    # signed volume, as Brainstorm stores them; otherwise every gain changes sign.
    if signed_volume(V, F) > 0:
        F = F[:, [0, 2, 1]]
        print(f"  {name}: faces re-ordered for OpenMEEG")
    return F


def is_closed(F):
    edges = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    return np.all(counts == 2)


def points_inside(P, V, F):
    """Winding-number test: True where P lies inside the closed mesh (V, F)."""
    inside = np.zeros(len(P), dtype=bool)
    A, B, C = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    for i0 in range(0, len(P), 500):
        p = P[i0:i0 + 500, None, :]
        a, b, c = A[None] - p, B[None] - p, C[None] - p
        la, lb, lc = (np.linalg.norm(x, axis=2) for x in (a, b, c))
        num = np.einsum("pfi,pfi->pf", a, np.cross(b, c))
        den = (la * lb * lc + np.einsum("pfi,pfi->pf", a, b) * lc
               + np.einsum("pfi,pfi->pf", b, c) * la + np.einsum("pfi,pfi->pf", c, a) * lb)
        inside[i0:i0 + 500] = np.abs(2 * np.arctan2(num, den).sum(axis=1)) > 2 * np.pi
    return inside


def read_channel_file(f):
    """Labels and positions (m) from a Brainstorm channel file or a label,x,y,z CSV."""
    if f.lower().endswith(".mat"):
        m = sio.loadmat(f, squeeze_me=True, struct_as_record=False)
        chans = [c for c in np.atleast_1d(m["Channel"])
                 if np.size(c.Loc) >= 3 and str(getattr(c, "Type", "EEG")).upper() == "EEG"]
        labels = [str(c.Name) for c in chans]
        # Loc is 3 x K (one column per point; the first is the electrode): read
        # column-wise (row-wise reading mixed the coordinates of K > 1)
        pos = np.array([np.asarray(c.Loc, dtype=float).reshape(3, -1)[:, 0] for c in chans])
    else:
        with open(f, newline="") as fh:
            rows = [{k.strip().lower(): v for k, v in r.items()} for r in csv.DictReader(fh)]
        labels = [r["label"].strip() for r in rows]
        pos = np.array([[float(r["x"]), float(r["y"]), float(r["z"])] for r in rows])
    if len(set(l.lower() for l in labels)) != len(labels):
        raise RuntimeError(f"duplicate labels in {f}")
    return labels, pos


def surface_points(V, F, n=6):
    """Points covering the triangles of a surface (barycentric grid of n steps per edge):
    distances to them are distances to the surface within about a millimetre, not to
    its vertices (about 1.5 cm apart on tess_head)."""
    pts = [V]
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            if (i, j, k).count(0) >= 2:        # the vertices themselves
                continue
            pts.append((i * V[F[:, 0]] + j * V[F[:, 1]] + k * V[F[:, 2]]) / n)
    return np.vstack(pts)


def check_positions(labels, elec, scalp_vertices, force=False, scalp_faces=None):
    """Stop when the electrode positions do not fit the template head: units,
    distance to the scalp, azimuth of the landmarks (x anterior, y left).
    With scalp_faces, the distance is to the scalp surface; without, to its
    vertices (or the points given)."""
    problems = []
    if np.abs(elec).max() > 1.0:
        problems.append(f"coordinates up to {np.abs(elec).max():.3g}: expected metres (mm or cm given?)")
    scalp = scalp_vertices if scalp_faces is None else surface_points(scalp_vertices, scalp_faces)
    d_el, _ = cKDTree(scalp).query(elec)
    print(f"{len(labels)} electrodes: distance to the scalp median {np.median(d_el) * 1e3:.1f} mm, "
          f"max {d_el.max() * 1e3:.1f} mm")
    if np.median(d_el) > 0.01:
        problems.append(f"electrodes {np.median(d_el) * 1e3:.0f} mm (median) from the template scalp")
    lower = [l.lower() for l in labels]
    az = np.degrees(np.arctan2(elec[:, 1], elec[:, 0]))
    landmarks = [("Fpz", 0), ("T7", 90), ("Oz", 180), ("T8", -90)]
    for l, expect in landmarks:
        if l.lower() in lower:
            a = az[lower.index(l.lower())]
            off = (a - expect + 180) % 360 - 180
            print(f"  {l}: azimuth {a:.0f} deg (expected about {expect})")
            if abs(off) > 30:
                problems.append(f"{l} at azimuth {a:.0f} deg, expected about {expect} (rotated cap?)")
    if not any(l.lower() in lower for l, _ in landmarks):
        print("  no Fpz, T7, Oz or T8 label: the orientation of the cap cannot be checked from its labels "
              "(check that x points to the nose and y to the left ear)")
    for p in problems:
        if force:
            print("  WARNING: " + p)
    if problems and not force:
        raise RuntimeError("electrode positions do not fit the template head: " + "; ".join(problems)
                           + " (check the coordinates, or --force)")
    return problems


def check_electrodes(labels, elec, max_dist=0.02):
    """Electrodes of the physics check: Cz, Fpz, Oz, T7, T8 by label, or, for a cap
    without these labels (e.g. ANT's equidistant Z1, LL1, ...), the electrode nearest
    to each one's position on the template (ASA 10-05 cap, within max_dist m).
    Returns (index, name) pairs."""
    lower = [l.lower() for l in labels]
    asa_labels, asa_pos = read_channel_file(os.path.join(ANAT, "channel_ASA_10-05_343_bs.mat"))
    asa = {l.lower(): p for l, p in zip(asa_labels, asa_pos)}
    tree = cKDTree(elec)
    out = []
    for lab in ["Cz", "Fpz", "Oz", "T7", "T8"]:
        if lab.lower() in lower:
            out.append((lower.index(lab.lower()), lab))
            continue
        d, ie = tree.query(asa[lab.lower()])
        if d <= max_dist:
            out.append((int(ie), f"{labels[ie]} (at {lab}, {d * 1e3:.0f} mm)"))
    return out


def as_brainstorm_writes(V):
    """Vertices (m) as Brainstorm writes them for OpenMEEG (out_tess_tri: millimetres, %g)."""
    return np.vectorize(lambda x: float("%g" % x))(V * 1000.0) / 1000.0


def gain(surfaces, dipoles, elec):
    """OpenMEEG EEG gain (nChannels x 3*nDipoles, x/y/z interleaved) of the nested
    surfaces [(V, F) inner skull, outer skull, scalp]."""
    geom = om.make_nested_geometry(surfaces, list(CONDUCTIVITY))
    if not geom.is_nested():
        raise RuntimeError("geometry not nested")
    nV = len(dipoles)
    dip = om.Matrix(np.asfortranarray(np.c_[np.repeat(dipoles, 3, axis=0), np.tile(np.eye(3), (nV, 1))]))
    sensors = om.Sensors(om.Matrix(np.asfortranarray(elec)), geom)
    t0 = datetime.datetime.now()
    hm = om.HeadMat(geom)
    hm.invert()
    G = np.array(om.GainEEG(hm, om.DipSourceMat(geom, dip, "Brain"), om.Head2EEGMat(geom, sensors)).array())
    print(f"Gain {G.shape} in {(datetime.datetime.now() - t0).seconds} s", flush=True)
    return G


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--channel-file", required=True, help="Brainstorm channel .mat or label,x,y,z CSV (m)")
    ap.add_argument("--out", required=True, help="output .mat")
    ap.add_argument("--free", action="store_true", help="also save G_free (x, y, z dipoles)")
    ap.add_argument("--force", action="store_true", help="position problems are warnings, not errors")
    ap.add_argument("--anatomy", default=ANAT,
                    help="folder of tess_head.mat, tess_outerskull.mat, tess_innerskull.mat (default: inputs/, ICBM152)")
    ap.add_argument("--cortex", default=None,
                    help="cortex surface of the sources (default: tess_cortex_pial_low.mat of inputs/, "
                         "or tess_cortex.mat of --anatomy)")
    args = ap.parse_args()
    anat = os.path.abspath(args.anatomy)
    cortex_file = args.cortex
    if cortex_file is None:
        same = os.path.normcase(anat) == os.path.normcase(os.path.abspath(ANAT))
        cortex_file = os.path.join(anat, "tess_cortex_pial_low.mat" if same else "tess_cortex.mat")

    Vh, Fh, _ = load_surface("tess_head.mat", anat)
    Vo, Fo, _ = load_surface("tess_outerskull.mat", anat)
    Vi, Fi, _ = load_surface("tess_innerskull.mat", anat)
    Vc, _, cortex = load_surface(os.path.basename(cortex_file), os.path.dirname(os.path.abspath(cortex_file)))
    Nc = np.asarray(cortex["VertNormals"], dtype=float)
    Nc /= np.linalg.norm(Nc, axis=1, keepdims=True)
    Fh, Fo, Fi = (orient_for_openmeeg(V, F, n) for V, F, n in
                  [(Vh, Fh, "head"), (Vo, Fo, "outer skull"), (Vi, Fi, "inner skull")])
    for n, F in [("head", Fh), ("outer skull", Fo), ("inner skull", Fi)]:
        if not is_closed(F):
            raise RuntimeError(f"{n} surface is not closed")
    assert points_inside(Vo, Vh, Fh).all(), "outer skull not inside scalp"
    assert points_inside(Vi, Vo, Fo).all(), "inner skull not inside outer skull"
    # A few pial vertices of the template poke through the inner skull: move each
    # 2 mm inside, towards the skull centroid, so that every dipole is in the brain.
    moved = np.flatnonzero(~points_inside(Vc, Vi, Fi))
    Vc_bem = Vc.copy()
    if moved.size:
        _, near = cKDTree(Vi).query(Vc[moved])
        inward = Vi.mean(axis=0) - Vi[near]
        inward /= np.linalg.norm(inward, axis=1, keepdims=True)
        Vc_bem[moved] = Vi[near] + 0.002 * inward
        print(f"  {moved.size} cortical vertices moved inside the inner skull: {moved + 1}")

    labels, elec = read_channel_file(args.channel_file)
    check_positions(labels, elec, Vh, args.force, scalp_faces=Fh)

    nV = len(Vc)
    surfaces = [(Vi, Fi), (Vo, Fo), (Vh, Fh)]
    G = gain(surfaces, Vc_bem, elec)
    precision = "full (double)"
    if not np.isfinite(G).all():
        # OpenMEEG's analytic integrals are not finite for some exactly degenerate
        # configurations of a mesh (e.g. Brainstorm's 1922-vertex ICBM152 skulls:
        # a face with an exactly coplanar neighbour). Brainstorm passes the
        # surfaces to OpenMEEG in millimetres with 6 significant digits
        # (out_tess_tri, %g), which removes them: the same here.
        print(f"  the gain has {np.sum(~np.isfinite(G))} non-finite values with the surfaces at full precision: "
              "again with their vertices as Brainstorm writes them for OpenMEEG (mm, %g)")
        surfaces = [(as_brainstorm_writes(V), F) for V, F in surfaces]
        G = gain(surfaces, Vc_bem, elec)
        precision = "as Brainstorm writes them for OpenMEEG (mm, %g): full precision gave non-finite values"
        if not np.isfinite(G).all():
            raise RuntimeError("the gain has non-finite values (OpenMEEG): check the BEM surfaces in Brainstorm")
    G_normal = np.einsum("cvk,vk->cv", G.reshape(len(elec), nV, 3), Nc)

    # physics check: a superficial dipole pointing at an electrode gives a positive,
    # near-maximal (average-referenced) potential there
    Gn = G_normal - G_normal.mean(axis=0, keepdims=True)
    n_checked = 0
    for ie, lab in check_electrodes(labels, elec):
        dist = np.linalg.norm(Vc - elec[ie], axis=1)
        u = (elec[ie] - Vc) / dist[:, None]
        cand = np.flatnonzero((dist < 0.035) & (np.einsum("ij,ij->i", Nc, u) > 0.9))
        if cand.size == 0:
            continue
        v = cand[np.argmin(dist[cand])]
        best = int(np.argmax(Gn[:, v]))
        # maximal there or at a neighbour (one step of the 10-05 system, 15-20 mm): the nearest
        # source facing the electrode can be off-centre (no vertices on the midline of an
        # individual cortex); a wrong sign, a rotated or misplaced cap fail
        ok = Gn[ie, v] > 0 and (best == ie or np.linalg.norm(elec[best] - elec[ie]) < 0.020)
        print(f"  check {lab}: V = {Gn[ie, v]:+.3g}, maximal at {labels[best]} -> {'ok' if ok else 'FAILED'}")
        n_checked += 1
        if not ok:
            raise RuntimeError(f"sign/localisation check failed at {lab}")
    if n_checked == 0:
        print("  WARNING: no electrode near Cz, Fpz, Oz, T7 or T8: sign/localisation not checked")

    out = {
        "G_normal": G_normal.astype(np.float32),
        "labels": np.array(labels, dtype=object),
        "elec_pos": elec, "src_pos": Vc, "src_nrm": Nc, "src_moved_idx": moved + 1,
        "conductivity": np.array(CONDUCTIVITY), "unit_gain": "V/(A*m)",
        "reference": "none (average reference applied by the pipeline)",
        "cap_source": os.path.basename(args.channel_file),
        "anatomy": f"{anat}: tess_head ({len(Vh)} V), tess_outerskull ({len(Vo)} V), tess_innerskull ({len(Vi)} V), "
                   f"{os.path.basename(cortex_file)} ({len(Vc)} V)",
        "bem_vertices": precision,
        "software": f"OpenMEEG {om.__version__} (Python), numpy {np.__version__}",
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    if args.free:
        out["G_free"] = G
    sio.savemat(args.out, out, do_compression=True)
    print("Saved", args.out, f"({os.path.getsize(args.out) / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
