"""Tests of the checks of tools/build_leadfield_openmeeg.py (V3-10).

Run with the Python that has OpenMEEG (README, section 1: What you need), from
v_dynamic/:  python -m pytest tests/python   or   python tests/python/test_build_leadfield.py
"""
import os
import sys
import tempfile

import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import build_leadfield_openmeeg as bl  # noqa: E402


def sphere_cap():
    """Electrodes on a 9 cm sphere and the sphere's points as the 'scalp' (x anterior, y left)."""
    labels = ["Fpz", "T7", "Oz", "T8", "Cz"]
    d = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1]], dtype=float)
    u = np.random.default_rng(0).normal(size=(4000, 3))
    scalp = 0.09 * u / np.linalg.norm(u, axis=1, keepdims=True)
    return labels, 0.09 * d, scalp


def raises(f):
    try:
        f()
    except RuntimeError:
        return True
    return False


def test_position_checks_are_errors():
    labels, elec, scalp = sphere_cap()
    assert bl.check_positions(labels, elec, scalp) == []
    assert raises(lambda: bl.check_positions(labels, elec * 1000, scalp * 1000))      # mm: the units check
    rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)                   # 90 degrees about z
    assert raises(lambda: bl.check_positions(labels, elec @ rot.T, scalp))
    assert raises(lambda: bl.check_positions(labels, elec * 1.3, scalp))              # 2.7 cm off the scalp
    assert len(bl.check_positions(labels, elec @ rot.T, scalp, force=True)) == 4      # --force: warnings


def test_distance_to_the_surface_not_its_vertices():
    # one large triangle; an electrode on its centre is ~0 from the surface, ~5 cm from its vertices
    V = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]])
    F = np.array([[0, 1, 2]])
    e = np.array([[0.1 / 3, 0.1 / 3, 0.0]])
    assert raises(lambda: bl.check_positions(["E1"], e, V))                 # vertices: 47 mm, an error
    assert bl.check_positions(["E1"], e, V, scalp_faces=F) == []            # surface: 0 mm


def test_check_electrodes_by_label_or_position():
    asa_labels, asa_pos = bl.read_channel_file(os.path.join(bl.ANAT, "channel_ASA_10-05_343_bs.mat"))
    pos = {l: p for l, p in zip(asa_labels, asa_pos)}
    names = ["Cz", "Fpz", "Oz", "T7", "T8"]
    # 10-20 labels: by label
    got = bl.check_electrodes(names, np.array([pos[n] for n in names]))
    assert [i for i, _ in got] == [0, 1, 2, 3, 4] and [n for _, n in got] == names
    # other labels, 5 mm from the landmarks: the nearest electrode, named with the landmark
    shift = np.array([0.005, 0.0, 0.0])
    other = ["Z5", "Z1", "Z11", "LL9", "RR9"]
    got = bl.check_electrodes(other, np.array([pos[n] + shift for n in names]))
    assert [i for i, _ in got] == [0, 1, 2, 3, 4]
    assert got[0][1].startswith("Z5 (at Cz")
    # nothing within 2 cm: no check
    assert bl.check_electrodes(other, np.array([pos[n] * 0.5 for n in names])) == []


def test_loc_read_column_wise():
    # a channel whose Loc is 3 x 2 (electrode and a second point): the first column
    Loc = np.array([[0.01, 0.07], [0.02, 0.08], [0.03, 0.09]])
    ch = [{"Name": "A", "Type": "EEG", "Loc": Loc}, {"Name": "B", "Type": "EEG", "Loc": np.array([0.04, 0.05, 0.06])}]
    f = os.path.join(tempfile.mkdtemp(), "channel_test.mat")
    sio.savemat(f, {"Channel": np.array([(c["Name"], c["Type"], c["Loc"]) for c in ch],
                                        dtype=[("Name", object), ("Type", object), ("Loc", object)])})
    labels, pos = bl.read_channel_file(f)
    assert labels == ["A", "B"]
    assert np.allclose(pos, [[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]])


def test_vertices_as_brainstorm_writes_them():
    # Brainstorm's out_tess_tri: millimetres, %g (6 significant digits)
    V = np.array([[0.12934123456, -0.0015896643744781613, 0.0291], [0.0, 1e-7, -0.08999999999]])
    W = bl.as_brainstorm_writes(V)
    assert np.array_equal(W, np.array([[129.341, -1.58966, 29.1], [0.0, 1e-4, -90.0]]) / 1000)
    assert np.abs(W - V).max() < 1e-6


if __name__ == "__main__":
    test_position_checks_are_errors()
    test_distance_to_the_surface_not_its_vertices()
    test_check_electrodes_by_label_or_position()
    test_loc_read_column_wise()
    test_vertices_as_brainstorm_writes_them()
    print("5 tests passed")
